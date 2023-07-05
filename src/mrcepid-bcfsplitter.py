#!/usr/bin/env python
#
# mrcepid-bcfsplitter
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/
import os
import csv
import logging
import dxpy

from pathlib import Path
from typing import List, Dict

from general_utilities.association_resources import run_cmd, generate_linked_dx_file, download_dxfile_by_name
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()


def ingest_docker_file() -> None:
    """Download the default Docker image so that we can run tools not on the DNANexus platform.

    :return: None
    """
    cmd = "docker pull egardner413/mrcepid-burdentesting:latest"
    run_cmd(cmd, is_docker=False)


def download_vcf(input_vcf: str) -> str:
    """Helper function that downloads a VCF, and it's index, then returns the filename prefix for that vcf

    :param input_vcf: An input_vcf file in the form of a DNANexus file ID (file-123456...)
    :return: A prefix of the downloaded VCF (missing .vcf.gz)
    """

    # Unsure if this is safe, but I think is the easiest way to avoid hardcoding a given project ID
    project_ID = dxpy.PROJECT_CONTEXT_ID

    # Set a DX file handler for the VCF file & index chunk
    # Both the DXFile class and download_dxfile method require a project_ID as they directly access UKBB bulk data.
    # The method I have used _hopefully_ avoids hardcoding...
    vcf = dxpy.DXFile(input_vcf, project=project_ID)

    # Now we need to find the corresponding tbi index using dxpy search functions
    idx_folder = vcf.describe()['folder']
    idx_name = vcf.describe()['name'] + '.tbi'
    idx_object = dxpy.find_one_data_object(more_ok=False, classname='file', project=project_ID, folder=idx_folder,
                                           name=idx_name, name_mode='exact')
    vcfidx = dxpy.DXFile(dxid=idx_object['id'], project=idx_object['project'])

    # And download both...
    download_dxfile_by_name(vcf, print_status=False)
    download_dxfile_by_name(vcfidx, print_status=False)

    # Set a prefix name for all files so that we can output a standard-named file:
    vcfprefix = vcf.describe()['name'].replace('.vcf.gz', '')

    return vcfprefix


def generate_variant_list(vcfprefix: str) -> None:
    """Generates a txt file of all variants in the provided BCF file

    This method is a wrapper around bcftools query to print a .tsv file with columns of CHROM and POS. This file can
    then be used to generate split lists of files to extract variants into seperate chunks.

    :param vcfprefix: A prefix of the downloaded VCF (missing .vcf.gz)
    :return: None
    """

    cmd = f'bcftools query -f "%CHROM\\t%POS\\n" -o /test/{vcfprefix}.bcf_sites.txt /test/{vcfprefix}.vcf.gz'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')


def split_sites(vcfprefix: str) -> List[str]:
    """Split site lists from a single VCF into multiple site lists to facilitate variant extraction

    This method is essentially is a re-code of the *NIX split command with a few special tweaks to:
    1. Ensure that BCF files smaller than 2,500 variants do not get created
    2. Make sure that variants do not get duplicated at the beginning/end of VCFs due to poor handling of MNVs

    :param vcfprefix: A prefix of the downloaded VCF (missing .vcf.gz)
    :return: A List of individual file chunk names of split variant lists
    """

    # Get number of variants in the sites file – this defines how many files we will write
    with Path(f'{vcfprefix}.bcf_sites.txt').open('r') as sites_reader:
        n_lines = 0
        for _ in sites_reader:
            n_lines += 0

        if n_lines == 0:
            raise ValueError(f'Could not determine file length of {vcfprefix}.bcf_sites.txt.')

        if n_lines % 5000 < 2500:
            append_last = True
        else:
            append_last = False

    # Then do second iteration to write the actual index files:
    with Path(f'{vcfprefix}.bcf_sites.txt').open('r') as sites_reader:
        # These variables are to control iteration parameters
        current_index_num = 1
        current_end = 0
        write_next_variant = True

        # Store the actual names of each chunk, so we can do the actual bcftools extraction later
        file_chunk_names = [f'{vcfprefix}_chunk{current_index_num}']

        # Collection of files that store different information for i/o:
        # 'sites' is the reader of ALL variants in the BCF
        # 'current_index' stores the actual variants to be extracted for each chunk
        sites = csv.DictReader(sites_reader,
                               delimiter='\t',
                               fieldnames=['chrom', 'pos'],
                               quoting=csv.QUOTE_NONE)
        current_index = Path(f'{vcfprefix}_chunk{current_index_num}').open('w')

        # Do the iteration itself
        for site in sites:
            if sites.line_num == 1:
                pass
            elif (sites.line_num - 1) % 5000 == 0:
                if (n_lines - sites.line_num) < 2500 and append_last is True:
                    logging.info("Appending remaining lines to end of last chunk...")
                else:
                    current_index.close()
                    current_index_num += 1
                    current_index = Path(f'{vcfprefix}_chunk{current_index_num}').open('w')
                    file_chunk_names.append(f'{vcfprefix}_chunk{current_index_num}')
                    if current_end == site['pos']:
                        write_next_variant = False
            # Allows me to control writing of variants
            if write_next_variant:
                current_index.write(f"{site['chrom']}\t{site['pos']}\n")
                current_end = site['pos']
            else:
                write_next_variant = True
        current_index.close()

    return file_chunk_names


def split_bcfs(vcfprefix: str, file_chunk_names: List[str]) -> List[dxpy.DXFile]:
    """A wrapper around bcftools view that generates smaller chunks from the original VCF file

    :param vcfprefix: A prefix of the downloaded VCF (missing .vcf.gz)
    :param file_chunk_names: A List of individual file chunk names of split variant lists
    :return: A List of DNANexus DXFile objects that have been uploaded to the DNANexus platform and are ready for
        return to the requesting project.
    """

    current_chunk = 1
    bcf_files = []
    for file in file_chunk_names:
        cmd = f'bcftools view --threads 4 -T /test/{file} ' \
              f'-Ob -o /test/{vcfprefix}_chunk{current_chunk}.bcf ' \
              f'/test/{vcfprefix}.vcf.gz'
        run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')
        bcf_files.append(generate_linked_dx_file(f'{vcfprefix}_chunk{current_chunk}.bcf'))
        current_chunk += 1
    Path(f'{vcfprefix}.vcf.gz').unlink()
    return bcf_files


def process_vcf(input_vcf: str) -> List[dxpy.DXFile]:
    """Helper function that enables multithreading in this applet

    This method calls the individual functions of the methods in this applet for a single VCF. This method just allows
    for calling each vcf as an individual thread within an instance by general_utilities.job_management.thread_utility
    and passing it through the splitting protocol outlined in the main method and README.

    :param input_vcf: A DNANexus file-id (file-12345...) in string format to process
    :return: A List of DNANexus DXFile objects that have been uploaded to the DNANexus platform and are ready for
        return to the requesting project created by the :func:`split_bcfs()` method.
    """

    # Download the VCF
    vcfprefix = download_vcf(input_vcf)

    # Do actual splitting of the target VCF:
    # 1. generate a list of all variants in the file
    generate_variant_list(vcfprefix)

    # 2. Generate reasonable sized (5k) lists of variants:
    file_chunk_names = split_sites(vcfprefix)

    # 3. And actually split the files into chunks:
    bcf_files = split_bcfs(vcfprefix, file_chunk_names)

    return bcf_files


@dxpy.entry_point('main')
def main(input_vcfs: dict, testing_script: dict, testing_directory: str) -> dict:
    """This is the :func:`main()` method for all apps/applets required by DNANexus.

    This applet is a fairly simple wrapper around bcftools. For each vcf.gz file provided to `input_vcfs` it:

    1. Prints a list of all sites found in the vcf via bcftools query
    2. Split sites into equal-sized chunks, if possible
    3. Extract these sites into individual .bcf files from the original vcf.gz using bcftools view -T

    :param input_vcfs: A DNANexus file-id (file-12345...) pointing to a list-file containing DNANexus file-ids
        (file-12345...) to split into smaller chunks.
    :param testing_script: Script compatible with pytest. If not null, invoke the bcfsplitter testing suite
        via :func:`test`.
    :param testing_directory: Directory containing test files if in testing mode.
    :return: A dictionary of outputs
    """

    if testing_script:
        LOGGER.info('Testing mode activated...')
        if testing_directory is None:
            raise ValueError(f'Testing mode invoked but -itesting_directory not provided!')

        output = test(input_vcfs, testing_script, testing_directory)
    else:
        # Ingest Docker file:
        ingest_docker_file()

        # Run through each VCF file provided and perform filtering.
        # input_vcfs is simple a file list of DNANexus file hashes that I dereference below
        input_vcfs = dxpy.DXFile(input_vcfs['$dnanexus_link'])
        dxpy.download_dxfile(input_vcfs.get_id(), 'vcf_list.txt')  # Actually download the file

        # Use thread utility to multi-thread this process:
        thread_utility = ThreadUtility(thread_factor=4,
                                       error_message='A splitting thread failed',
                                       incrementor=5)

        # And launch individual jobs
        with Path('vcf_list.txt').open('r') as input_vcf_reader:
            for line in input_vcf_reader:
                input_vcf = line.rstrip()
                thread_utility.launch_job(process_vcf,
                                          input_vcf=input_vcf)

        bcf_files = []
        for result in thread_utility:
            bcf_files.extend(result)

        # Set output
        output = {"output_vcfs": [dxpy.dxlink(item) for item in bcf_files]}

    return output


def test(input_vcfs: dict, testing_script: dict, testing_directory: str) -> Dict:
    """Run the bcfsplitter testing suite.

    This method is invisible to the applet and can only be accessed by using API calls via dxpy.DXApplet() on
    a local machine. See the resources in the `./test/` folder for more information on running tests.

    :param input_vcfs: A DNANexus file-id (file-12345...) pointing to a list-file containing DNANexus file-ids
        (file-12345...) to split into smaller chunks.
    :param testing_script: Script compatible with pytest. If not null, invoke the bcfsplitter testing suite
        via :func:`test`.
    :param testing_directory: Directory containing test files if in testing mode.
    :return:
    """

    LOGGER.info('Launching mrcepid-filterbcf with the testing suite')
    dxpy.download_dxfile(dxid=testing_script['$dnanexus_link'], filename='test.py')

    # I then set an environment variable that tells pytest where the testing directory is
    os.environ['CI'] = '500'  # Make sure pytest logs aren't truncated
    os.environ['TEST_DIR'] = testing_directory
    LOGGER.info(f'TEST_DIR environment variable set: {os.getenv("TEST_DIR")}')
    os.environ['INPUT_VCFS'] = input_vcfs['$dnanexus_link']
    LOGGER.info(f'INPUT_VCFS environment variable set: {os.getenv("BGEN_INDEX")}')

    out_log = Path(f'pytest.{output_prefix}.log')
    try:
        run_cmd('pytest test.py', is_docker=False, stdout_file=out_log)
    except RuntimeError:
        pass

    return {}


dxpy.run()
