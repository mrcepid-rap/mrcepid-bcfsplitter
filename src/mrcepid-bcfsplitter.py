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


def ingest_human_reference(human_reference: dict, human_reference_index: dict) -> None:
    """Download human reference files – default dxIDs are the location of the GRCh38 reference file on AWS London

    :param human_reference: DXLink to the human reference file (must be a .fa.gz)
    :param human_reference_index: DXLink to the human reference file index
    """

    dxpy.download_dxfile(dxpy.DXFile(human_reference).get_id(), "reference.fasta.gz")
    dxpy.download_dxfile(dxpy.DXFile(human_reference_index).get_id(), "reference.fasta.fai")
    cmd = "gunzip reference.fasta.gz"  # Better to unzip the reference for most commands for some reason...
    CMD_EXEC.run_cmd(cmd)


def replace_multi_suffix(original_path: Path, new_suffix: str) -> Path:
    """A helper function to replace a path on a file with multiple suffixes (e.g., .tsv.gz)

    This function just loops through the path and recursively removes the string after '.'. Once there are no more
    full stops it then adds the requested :param: new_suffix.

    :param original_path: The original filepath
    :param new_suffix: The new suffix to add
    :return: A Pathlike to the new file
    """

    while original_path.suffix:
        original_path = original_path.with_suffix('')

    return original_path.with_suffix(new_suffix)


def generate_site_tsv(vcf_file: Path, sites_suffix: str) -> Path:
    """This method is a helper to generate a sites file for a given input VCF file

    A simple wrapper around `bcftools query` to generate a .tsv file with columns of CHROM, POS, REF, ALT, AC. This file
    will not contain a header, and thus needs to be considered when reading the file with a :func:`csv.DictReader`.

    :param vcf_file: A Pathlike representation of a VCF or BCF file (can also be gzipped).
    :param sites_suffix: The suffix to append to the sites file.
    :return: A Pathlike to the generated sites file
    """

    sites_file = replace_multi_suffix(vcf_file, sites_suffix)

    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\n" ' \
          f'-o /test/{sites_file} /test/{vcf_file}'
    CMD_EXEC.run_cmd_on_docker(cmd)

    return sites_file


def download_vcf(input_vcf: str) -> Tuple[Path, int]:
    """Helper function that downloads a VCF and it's index, and then calculates the filesize for record-keeping purposes.

    :param input_vcf: An input_vcf file in the form of a DNANexus file ID (file-123456...)
    :return: The downloaded VCF file and the size (in bytes) of the VCF file
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

    # Get the size of the VCF file for logging purposes:
    vcf_size = vcfpath.stat().st_size

    return vcfpath, vcf_size


def count_variant_list_and_filter(vcf_file: Path, alt_allele_threshold: int) -> Tuple[Path, List[Dict], int, int, int]:
    """Counts the total number of variants in a given VCF/BCF and filters out sites with excessive alternate alleles

    This method 1st uses the :func:`generate_site_tsv` helper function to query the bcf and print a .tsv file.

    This method then filters sites with excessive alternate alleles. This is because to load them into memory would
    result in excessive memory use and cause the machine running it to crash. This filtered set of variants is then
    written to a sites file to be used for filtering of sites. Filtered sites are retained as a dictionary and
    returned, so they can be written to a file that documents filtered sites. Alleles with a minor allele frequency >
    0.1% are specifically captured in order to know blind spots during subsequent analysis.

    :param vcf_file: A vcf file to generate a sites list for
    :param alt_allele_threshold: Number of alternate alleles to allow in a variant before it is excluded
    :return: A tuple of the generated sites file as a Pathlike, a dict of failed sites, the number of variant rows, the
        number of pass variant rows, and the number of alternate alleles
    """

    sites_file = generate_site_tsv(vcf_file, '.sites.unfiltered.txt')
    filtered_sites_file = replace_multi_suffix(vcf_file, '.sites.txt')

    with sites_file.open('r') as sites_reader, \
            filtered_sites_file.open('w') as filtered_sites_writer:

        filtered_sites_csv = csv.DictWriter(filtered_sites_writer, delimiter='\t', extrasaction='ignore',
                                            fieldnames=['CHROM', 'POS', 'REF', 'ALT'])
        n_vcf_lines = 0
        n_final_lines = 0
        n_vcf_alternates = 0

        failed_sites = []

        # Note to future devs: DO NOT use a csv.DictReader for sites_reader! ALT alleles can be too long for the
        # standard buffer size and will cause a crash. This is why we use a standard reader here.
        for line in sites_reader:
            variant = dict(zip(['CHROM', 'POS', 'REF', 'ALT', 'AC'], line.rstrip().split('\t')))
            if variant['CHROM'] == 'CHROM':
                continue
            else:
                n_vcf_lines += 1
                alt_alleles = variant['ALT'].split(',')
                acs = list(map(int, variant['AC'].split(',')))
                n_var_alt = len(alt_alleles)

                # Check if the variant has too many alternate alleles
                if n_var_alt > alt_allele_threshold:

                    # Capture the site that failed
                    failed_site = {'chrom': variant['CHROM'], 'pos': variant['POS'], 'alts': [], 'acs': []}

                    # Capture high MAC alternate alleles at this site
                    for alt_num, alt in enumerate(alt_alleles):
                        if acs[alt_num] > 1000:  # Capture alts MAF ~> 0.1%
                            failed_site['alts'].append(alt)
                            failed_site['acs'].append(acs[alt_num])

                    failed_site['alts'] = ','.join(failed_site['alts'])
                    failed_site['acs'] = ','.join(map(str, failed_site['acs']))

                    failed_sites.append(failed_site)
                else:
                    filtered_sites_csv.writerow(variant)
                    n_final_lines += 1
                    n_vcf_alternates += n_var_alt

        return filtered_sites_file, failed_sites, n_vcf_lines, n_final_lines, n_vcf_alternates


def normalise_and_left_correct(bcf_file: Path, site_list: Path) -> Path:
    """A wrapper for BCFtools norm to left-normalise and split all variants

    Generate a normalised bcf file for all downstream processing using `bcftools norm`:
    -m : splits all multiallelics into separate records
    -f : provides a reference file so bcftools can left-normalise and check records against the reference genome
    -T : restricts the normalisation to a list of sites given by :param site_list:
    --old-rec-tag : sets a tag in the resulting bcf that contains the original record – used for IDing multi-allelics
        after splitting

    :param bcf_file: A list of bcffiles to normalise and left-correct
    :param site_list: A list of sites to restrict to and then normalise and left-correct
    :return: A Path object to the normalised and left-corrected BCF file
    """

    out_bcf = replace_multi_suffix(bcf_file, '.norm.bcf')
    cmd = f'bcftools norm --threads 8 -w 100 -Ob -m - -f /test/reference.fasta ' \
          f'-T /test/{site_list} ' \
          f'--old-rec-tag MA ' \
          f'-o /test/{out_bcf} /test/{bcf_file}'
    CMD_EXEC.run_cmd_on_docker(cmd)
    bcf_file.unlink()

    return out_bcf


def split_sites(filtered_sites_file: Path, n_lines: int, chunk_size: int) -> List[Path]:
    """Split site lists from a single VCF into multiple site lists to facilitate variant extraction

    This method is essentially is a re-code of the *NIX split command with a few special tweaks to:
    1. Ensure that BCF files smaller than chunk_size / 2 variants do not get created
    2. Make sure that variants do not get duplicated at the beginning/end of VCFs due to poor handling of multi-allelics

    :param filtered_sites_file: The filtered sites file
    :param n_lines: Number of sites in the provided BCF file
    :param chunk_size: Number of sites in the resulting output BCF file
    :return: A List of individual file chunk Pathlikes of split variant lists
    """

    # Determine if we have to append or create a new file based on line number cutoff
    if n_lines % chunk_size < (chunk_size / 2):
        append_last = True
    else:
        append_last = False

    # Then do second iteration to write the actual index files:
    with filtered_sites_file.open('r') as sites_reader:
        # These variables are to control iteration parameters
        current_index_num = 1
        current_end = 0
        write_next_variant = True

        # Store the actual names of each chunk, so we can do the actual bcftools extraction later
        current_index_path = replace_multi_suffix(filtered_sites_file,
                                                  f'.chunk{current_index_num}')
        file_chunk_paths = [current_index_path]

        # Collection of files that store different information for i/o:
        # 'sites' is the reader of ALL variants in the BCF
        # 'current_index' stores the actual variants to be extracted for each chunk
        sites = csv.DictReader(sites_reader,
                               delimiter='\t',
                               fieldnames=['chrom', 'pos', 'ref', 'alt', 'ac'],
                               quoting=csv.QUOTE_NONE)
        current_index_writer = current_index_path.open('w')

        # Do the iteration itself
        for site in sites:
            if sites.line_num == 1:
                pass
            elif (sites.line_num - 1) % chunk_size == 0:
                if (n_lines - sites.line_num) < (chunk_size / 2) and append_last is True:
                    logging.info("Appending remaining lines to end of last chunk...")
                else:
                    current_index_writer.close()
                    current_index_num += 1
                    current_index_path = replace_multi_suffix(filtered_sites_file,
                                                              f'.chunk{current_index_num}')
                    current_index_writer = current_index_path.open('w')
                    file_chunk_paths.append(current_index_path)

            # Don't write duplicate positions as BCFTools can only process on position rather than position / ref / alt
            if int(site['pos']) == current_end:
                write_next_variant = False

            # Allows me to control writing of variants
            if write_next_variant:
                current_index_writer.write(f"{site['chrom']}\t{site['pos']}\n")
                current_end = int(site['pos'])
            else:
                write_next_variant = True

        # Make sure we close the final writer
        current_index_writer.close()

    return file_chunk_paths


def split_bcfs(vcf_file: Path, file_chunk_names: List[Path]) -> List[dxpy.DXFile]:
    """A wrapper around bcftools view that generates smaller chunks from the original VCF file

    We also remove star '*' alleles at this point as we do not use them for any subsequent association testing.

    :param vcf_file: Path to the downloaded VCF to split into chunks
    :param file_chunk_names: A List of individual file chunk Paths of split variant lists
    :return: A List of individual file chunk Pathlikes of split BCF files
    """

    bcf_files = []
    for file_chunk in file_chunk_names:
        out_bcf = file_chunk.with_suffix(f'{file_chunk.suffix}.bcf')
        cmd = f'bcftools view --threads 8 -e "alt==\'*\'" -T /test/{file_chunk} ' \
              f'-Ob -o /test/{out_bcf} ' \
              f'/test/{vcf_file}'
        CMD_EXEC.run_cmd_on_docker(cmd)
        bcf_files.append(generate_linked_dx_file(out_bcf))

    vcf_file.unlink()

    return bcf_files


def process_vcf(input_vcf: str, chunk_size: int, alt_allele_threshold: int) -> Tuple[List[dxpy.DXFile], Dict, List[Dict]]:
    """Helper function that enables multithreading in this applet

    This method calls the individual functions of the methods in this applet for a single VCF. This method just allows
    for calling each vcf as an individual thread within an instance by general_utilities.job_management.thread_utility
    and passing it through the splitting protocol outlined in the main method and README.

    :param input_vcf: A DNANexus file-id (file-12345...) in string format to process
    :param chunk_size: The number of variants to include per-output BCF produced by this applet. Lines per-file cannot
        be smaller than [chunk_size] / 2.
    :param alt_allele_threshold: Number of alternate alleles to allow in a variant before it is excluded
    :return: A Tuple consisting of a List of DNANexus DXFile objects that have been uploaded to the DNANexus platform
        and are ready for return to the requesting project created by the :func:`split_bcfs()` method and a dict of
        information about the file that has been split for documentation purposes.
    """

    # 1. Download the VCF
    vcf_path, vcf_size = download_vcf(input_vcf)

    # 2. generate a list of all variants in the file, filtering for large alt size
    orig_sites, failed_sites, n_norm_lines, n_final_lines, n_norm_alts = count_variant_list_and_filter(vcf_path,
                                                                                                       alt_allele_threshold)

    # Collate information about this file
    log_info = {'vcf': vcf_path.name, 'dxid': input_vcf, 'n_sites': n_norm_lines, 'n_final_sites': n_final_lines, 'vcf_size': vcf_size}

    # Some WGS-based vcfs are meant to have 0 sites, and we want to capture that here, so we have a full accounting
    if n_norm_lines == 0:
        final_files = []  # Empty since no splitting will happen

    else:

        # 3. Normalise and left-correct all variants
        norm_bcf = normalise_and_left_correct(vcf_path, orig_sites)

        norm_sites = generate_site_tsv(norm_bcf, '.norm.sites.txt')

        # 4. Generate reasonable sized (param: chunk_size) lists of variants:
        file_chunk_paths = split_sites(norm_sites, n_norm_lines, chunk_size)

        # 5. And actually split the files into chunks:
        final_files = split_bcfs(norm_bcf, file_chunk_paths)

    return final_files, log_info, failed_sites


@dxpy.entry_point('main')
def main(input_vcfs: dict, chunk_size: int, alt_allele_threshold: int, output_name: str, human_reference: dict,
         human_reference_index: dict, testing_script: dict, testing_directory: str) -> dict:
    """This is the :func:`main()` method for all apps/applets required by DNANexus.

    This applet is a fairly simple wrapper around bcftools. For each vcf.gz file provided to `input_vcfs` it:

    1. Prints a list of all sites found in the vcf via bcftools query
    2. Split sites into equal-sized chunks, if possible
    3. Extract these sites into individual .bcf files from the original vcf.gz using bcftools view -T

    :param input_vcfs: A DNANexus file-id (file-12345...) pointing to a list-file containing DNANexus file-ids
        (file-12345...) to split into smaller chunks.
    :param chunk_size: The number of variants to include per-output BCF produced by this applet. Lines per-file cannot
        be smaller than [chunk_size] / 2.
    :param alt_allele_threshold: Number of alternate alleles to allow in a variant before it is excluded
    :param output_name: A prefix to use for output files
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
        thread_utility = ThreadUtility(thread_factor=8,
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
                                          chunk_size=chunk_size,
                                          alt_allele_threshold=alt_allele_threshold)

        bcf_files = []
        output_name = f'.{output_name}.' if output_name else '.'
        split_info_path = Path(f'vcf_info{output_name}tsv')
        skipped_sites_path = Path(f'skipped_sites{output_name}tsv')
        size_zero_bcf_count = 0
        with split_info_path.open('w') as split_info_file, \
                skipped_sites_path.open('w') as skipped_sites_file:
            split_info_csv = csv.DictWriter(split_info_file,
                                            fieldnames=['vcf', 'dxid', 'n_sites', 'n_final_sites', 'vcf_size'],
                                            delimiter='\t')
            split_info_csv.writeheader()

            skipped_sites_csv = csv.DictWriter(skipped_sites_file,
                                               fieldnames=['vcf', 'chrom', 'pos', 'alts', 'acs'],
                                               delimiter='\t')
            skipped_sites_csv.writeheader()

            for result in thread_utility:
                files, info, skipped_sites = result
                bcf_files.extend(files)
                if info['n_sites'] == 0:
                    size_zero_bcf_count += 1
                split_info_csv.writerow(info)
                for site in skipped_sites:
                    site['vcf'] = info['vcf']
                skipped_sites_csv.writerows(skipped_sites)

        LOGGER.info(f'Number of VCFs with 0 sites: {size_zero_bcf_count} '
                    f'({(size_zero_bcf_count / n_vcfs) * 100:0.2f}%)')

        # Set output
        output = {'output_vcfs': [dxpy.dxlink(item) for item in bcf_files],
                  'run_info': dxpy.dxlink(generate_linked_dx_file(split_info_path)),
                  'skipped_sites': dxpy.dxlink(generate_linked_dx_file(skipped_sites_path))}

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
