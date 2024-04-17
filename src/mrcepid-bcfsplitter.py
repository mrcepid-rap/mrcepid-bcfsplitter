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
import dxpy
import logging

from pathlib import Path
from typing import List, Dict, Tuple

from general_utilities.association_resources import generate_linked_dx_file, download_dxfile_by_name
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()
CMD_EXEC = build_default_command_executor()


def download_vcf(input_vcf: str) -> Tuple[str, int]:
    """Helper function that downloads a VCF, and it's index, then returns the filename prefix for that vcf

    :param input_vcf: An input_vcf file in the form of a DNANexus file ID (file-123456...)
    :return: A prefix of the downloaded VCF (missing .vcf.gz)
    """

    # Unsure if this is safe, but I think is the easiest way to avoid hardcoding a given project ID
    project_id = dxpy.PROJECT_CONTEXT_ID

    # Set a DX file handler for the VCF file & index chunk
    # Both the DXFile class and download_dxfile method require a project_id as they directly access UKBB bulk data.
    # The method I have used _hopefully_ avoids hardcoding...
    vcf = dxpy.DXFile(input_vcf, project=project_id)

    # Now we need to find the corresponding tbi index using dxpy search functions
    idx_folder = vcf.describe()['folder']
    idx_name = vcf.describe()['name'] + '.tbi'
    idx_object = dxpy.find_one_data_object(more_ok=False, classname='file', project=project_id, folder=idx_folder,
                                           name=idx_name, name_mode='exact')
    vcfidx = dxpy.DXFile(dxid=idx_object['id'], project=idx_object['project'])

    # And download both...
    vcfpath = download_dxfile_by_name(vcf, project_id=project_id, print_status=False)
    download_dxfile_by_name(vcfidx, project_id=project_id, print_status=False)

    # Set a prefix name for all files so that we can output a standard-named file:
    vcf_prefix = vcfpath.name.replace('.vcf.gz', '')
    vcf_size = vcfpath.stat().st_size

    return vcf_prefix, vcf_size


def ingest_human_reference(human_reference: dict, human_reference_index: dict) -> None:
    """Download human reference files – default dxIDs are the location of the GRCh38 reference file on AWS London

    :param human_reference: DXLink to the human reference file (must be a .fa.gz)
    :param human_reference_index: DXLink to the human reference file index
    """

    dxpy.download_dxfile(dxpy.DXFile(human_reference).get_id(), "reference.fasta.gz")
    dxpy.download_dxfile(dxpy.DXFile(human_reference_index).get_id(), "reference.fasta.fai")
    cmd = "gunzip reference.fasta.gz"  # Better to unzip the reference for most commands for some reason...
    CMD_EXEC.run_cmd(cmd)


def normalise_and_left_correct(vcf_prefix: str) -> None:
    """A wrapper for BCFtools norm to left-normalise and split all variants

    Generate a normalised bcf file for all downstream processing:
    -m : splits all multiallelics into separate records
    -f : provides a reference file so bcftools can left-normalise and check records against the reference genome
    --old-rec-tag : sets a tag in the resulting bcf that contains the original record – used for IDing multi-allelics
        after splitting
    """
    cmd = f'bcftools norm --threads 2 -w 500 -Ob -m - -f /test/reference.fasta ' \
          f'--old-rec-tag MA ' \
          f'-o /test/{vcf_prefix}.norm.bcf /test/{vcf_prefix}.vcf.gz'
    CMD_EXEC.run_cmd_on_docker(cmd)
    Path(f'{vcf_prefix}.vcf.gz').unlink()


def generate_and_count_variant_list(vcf_file: Path) -> Tuple[Path, int, int]:
    """Generates a txt file of all variants in the provided BCF file

    This method is a wrapper around bcftools query to print a .tsv file with columns of CHROM and POS. This file can
    then be used to generate split lists of files to extract variants into separate chunks.

    This method then iterates through the generated file

    :param vcf_file: A vcf file to generate a sites list for
    :return: A tuple of the generated sites file, the number of variant rows, and the number of alternate alleles
    """

    sites_file = Path(f'{vcf_file}.sites.txt')

    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" ' \
          f'-o /test/{sites_file} /test/{vcf_file}'
    CMD_EXEC.run_cmd_on_docker(cmd)

    with sites_file.open('r') as sites_reader:

        sites_csv = csv.DictReader(sites_reader, delimiter='\t', fieldnames=['CHROM', 'POS', 'REF', 'ALT'])
        n_vcf_lines = 0
        n_vcf_alternates = 0
        for variant in sites_csv:
            n_vcf_lines += 1
            n_var_alt = len(variant['ALT'].split(','))
            n_vcf_alternates += n_var_alt
            if n_var_alt >= 20:
                LOGGER.warning(f'Variant {variant["CHROM"]}:{variant["POS"]} from {vcf_file} has excessive number of '
                               f'alt alleles ({n_var_alt})')

        return sites_file, n_vcf_lines, n_vcf_alternates


def split_sites(vcf_prefix: str, n_lines: int, chunk_size: int) -> List[str]:
    """Split site lists from a single VCF into multiple site lists to facilitate variant extraction

    This method is essentially is a re-code of the *NIX split command with a few special tweaks to:
    1. Ensure that BCF files smaller than 2,500 variants do not get created
    2. Make sure that variants do not get duplicated at the beginning/end of VCFs due to poor handling of MNVs

    :param vcf_prefix: A prefix of the downloaded VCF (missing .vcf.gz)
    :param n_lines: Number of sites in the provided BCF file
    :param chunk_size: Number of sites in the resulting output BCF file
    :return: A List of individual file chunk names of split variant lists
    """

    # Determine if we have to append or create a new file based on line number cutoff
    if n_lines % chunk_size < (chunk_size / 2):
        append_last = True
    else:
        append_last = False

    # Then do second iteration to write the actual index files:
    with Path(f'{vcf_prefix}.norm.bcf.sites.txt').open('r') as sites_reader:
        # These variables are to control iteration parameters
        current_index_num = 1
        current_end = 0
        write_next_variant = True

        # Store the actual names of each chunk, so we can do the actual bcftools extraction later
        file_chunk_names = [f'{vcf_prefix}_chunk{current_index_num}']

        # Collection of files that store different information for i/o:
        # 'sites' is the reader of ALL variants in the BCF
        # 'current_index' stores the actual variants to be extracted for each chunk
        sites = csv.DictReader(sites_reader,
                               delimiter='\t',
                               fieldnames=['chrom', 'pos', 'ref', 'alt'],
                               quoting=csv.QUOTE_NONE)
        current_index = Path(f'{vcf_prefix}_chunk{current_index_num}').open('w')

        # Do the iteration itself
        for site in sites:
            if sites.line_num == 1:
                pass
            elif (sites.line_num - 1) % chunk_size == 0:
                if (n_lines - sites.line_num) < (chunk_size / 2) and append_last is True:
                    logging.info("Appending remaining lines to end of last chunk...")
                else:
                    current_index.close()
                    current_index_num += 1
                    current_index = Path(f'{vcf_prefix}_chunk{current_index_num}').open('w')
                    file_chunk_names.append(f'{vcf_prefix}_chunk{current_index_num}')

            # Don't write duplicate positions as BCFTools can only process on position rather than position / ref / alt
            if int(site['pos']) == current_end:
                write_next_variant = False

            # Allows me to control writing of variants
            if write_next_variant:
                current_index.write(f"{site['chrom']}\t{site['pos']}\n")
                current_end = int(site['pos'])
            else:
                write_next_variant = True
        current_index.close()

    return file_chunk_names


def split_bcfs(vcf_prefix: str, file_chunk_names: List[str]) -> List[dxpy.DXFile]:
    """A wrapper around bcftools view that generates smaller chunks from the original VCF file

    We also remove star '*' alleles at this point as we do not use them for any subsequent association testing.

    :param vcf_prefix: A prefix of the downloaded VCF (missing .vcf.gz)
    :param file_chunk_names: A List of individual file chunk names of split variant lists
    :return: A List of DNANexus DXFile objects that have been uploaded to the DNANexus platform and are ready for
        return to the requesting project.
    """

    current_chunk = 1
    bcf_files = []
    for file in file_chunk_names:
        cmd = f'bcftools view --threads 2 -e "alt==\'*\'" -T /test/{file} ' \
              f'-Ob -o /test/{vcf_prefix}_chunk{current_chunk}.bcf ' \
              f'/test/{vcf_prefix}.norm.bcf'
        CMD_EXEC.run_cmd_on_docker(cmd)
        bcf_files.append(generate_linked_dx_file(f'{vcf_prefix}_chunk{current_chunk}.bcf'))
        current_chunk += 1
    Path(f'{vcf_prefix}.norm.bcf').unlink()

    return bcf_files


def process_vcf(input_vcf: str, chunk_size: int) -> Tuple[List[dxpy.DXFile], Dict]:
    """Helper function that enables multithreading in this applet

    This method calls the individual functions of the methods in this applet for a single VCF. This method just allows
    for calling each vcf as an individual thread within an instance by general_utilities.job_management.thread_utility
    and passing it through the splitting protocol outlined in the main method and README.

    :param input_vcf: A DNANexus file-id (file-12345...) in string format to process
    :param chunk_size: The number of variants to include per-output BCF produced by this applet. Lines per-file cannot
        be smaller than [chunk_size] / 2.
    :return: A Tuple consisting of a List of DNANexus DXFile objects that have been uploaded to the DNANexus platform
        and are ready for return to the requesting project created by the :func:`split_bcfs()` method and a dict of
        information about the file that has been split for documentation purposes.
    """

    # Download the VCF
    vcf_prefix, vcf_size = download_vcf(input_vcf)

    # 2. Normalise and left-correct all variants
    normalise_and_left_correct(vcf_prefix)

    # Do actual splitting of the target VCF:
    # 3. generate a list of all variants in the file
    norm_sites, n_norm_lines, n_norm_alts = generate_and_count_variant_list(Path(f'{vcf_prefix}.norm.bcf'))

    # Collate information about this file
    log_info = {'vcf_prefix': vcf_prefix, 'n_sites': n_norm_lines, 'vcf_size': vcf_size}

    # Some WGS-based vcfs are meant to have 0 sites, and we want to capture that here, so we have a full accounting
    if n_norm_lines == 0:
        bcf_files = []  # Empty since no splitting will happen

    else:

        # 4. Generate reasonable sized (param: chunk_size) lists of variants:
        file_chunk_names = split_sites(vcf_prefix, n_norm_lines, chunk_size)

        # 5. And actually split the files into chunks:
        bcf_files = split_bcfs(vcf_prefix, file_chunk_names)

    return bcf_files, log_info


@dxpy.entry_point('main')
def main(input_vcfs: dict, chunk_size: int, human_reference: dict, human_reference_index: dict, testing_script: dict,
         testing_directory: str) -> dict:
    """This is the :func:`main()` method for all apps/applets required by DNANexus.

    This applet is a fairly simple wrapper around bcftools. For each vcf.gz file provided to `input_vcfs` it:

    1. Prints a list of all sites found in the vcf via bcftools query
    2. Split sites into equal-sized chunks, if possible
    3. Extract these sites into individual .bcf files from the original vcf.gz using bcftools view -T

    :param input_vcfs: A DNANexus file-id (file-12345...) pointing to a list-file containing DNANexus file-ids
        (file-12345...) to split into smaller chunks.
    :param chunk_size: The number of variants to include per-output BCF produced by this applet. Lines per-file cannot
        be smaller than [chunk_size] / 2.
    :param human_reference_index: Location of the human reference file in dxlink format
    :param human_reference: Location of the human reference file index in dxlink format
    :param testing_script: Script compatible with pytest. If not null, invoke the bcfsplitter testing suite
        via :func:`test`.
    :param testing_directory: Directory containing test files if in testing mode.
    :return: A dictionary of outputs as specified by dxapp.json
    """

    if testing_script:
        LOGGER.info('Testing mode activated...')
        if testing_directory is None:
            raise ValueError(f'Testing mode invoked but -itesting_directory not provided!')

        output = test(input_vcfs, testing_script, testing_directory)
    else:

        # Run through each VCF file provided and perform filtering.
        # input_vcfs is a file list of DNANexus file hashes that I dereference below
        input_vcfs = dxpy.DXFile(input_vcfs['$dnanexus_link'])
        dxpy.download_dxfile(input_vcfs.get_id(), 'vcf_list.txt')  # Actually download the file

        ingest_human_reference(human_reference, human_reference_index)

        # Use thread utility to multi-thread this process:
        thread_utility = ThreadUtility(thread_factor=2,
                                       error_message='A splitting thread failed',
                                       incrementor=5)

        # And launch individual jobs
        n_vcfs = 0
        with Path('vcf_list.txt').open('r') as input_vcf_reader:
            for line in input_vcf_reader:
                n_vcfs += 1
                input_vcf = line.rstrip()
                thread_utility.launch_job(process_vcf,
                                          input_vcf=input_vcf,
                                          chunk_size=chunk_size)

        bcf_files = []
        split_info_path = Path('vcf_info.tsv')
        size_zero_bcf_count = 0
        with split_info_path.open('w') as split_info_file:
            split_info_csv = csv.DictWriter(split_info_file,
                                            fieldnames=['vcf_prefix', 'n_sites', 'vcf_size'],
                                            delimiter='\t')
            split_info_csv.writeheader()

            for result in thread_utility:
                files, info = result
                bcf_files.extend(files)
                if info['n_sites'] == 0:
                    size_zero_bcf_count += 1
                split_info_csv.writerow(info)

        LOGGER.info(f'Number of VCFs with 0 sites: {size_zero_bcf_count} '
                    f'({(size_zero_bcf_count / n_vcfs)*100:0.2f}%)')
        if size_zero_bcf_count == 0:
            LOGGER.warning(f'All VCFs in this run were empty. This job will produce 0 output BCF files.')

        # Set output
        output = {'output_vcfs': [dxpy.dxlink(item) for item in bcf_files],
                  'run_info': dxpy.dxlink(generate_linked_dx_file(split_info_path))}

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

    out_log = Path(f'pytest.log')
    try:
        CMD_EXEC.run_cmd('pytest test.py', stdout_file=out_log)
    except RuntimeError:
        pass

    return {}


dxpy.run()
