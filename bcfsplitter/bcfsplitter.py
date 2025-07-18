#!/usr/bin/env python
#
# mrcepid-bcfsplitter
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import csv
import dxpy
import logging

from pathlib import Path
from typing import List, Dict, Tuple, Optional

from general_utilities.association_resources import replace_multi_suffix, find_index
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler, FileType
from general_utilities.import_utils.file_handlers.dnanexus_utilities import generate_linked_dx_file
from general_utilities.job_management.command_executor import build_default_command_executor, CommandExecutor
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()
CMD_EXEC = build_default_command_executor()


def ingest_human_reference(human_reference: dict, human_reference_index: dict, cmd_exec: CommandExecutor = CMD_EXEC) -> Path:
    """Download human reference files – default dxIDs are the location of the GRCh38 reference file on AWS London

    :param human_reference: DXLink to the human reference file (must be a .fa.gz)
    :param human_reference_index: DXLink to the human reference file index
    :param cmd_exec: An optional alternate CommandExecutor to use for executing system calls. Defaults to the default
        CMD_EXEC from this module.
    """

    reference = InputFileHandler(human_reference).get_file_handle()
    InputFileHandler(human_reference_index).get_file_handle()

    # if the file is not unzipped, unzip it
    if reference.suffix == '.gz':
        cmd = f"gunzip {reference}"  # Better to unzip the reference for most commands for some reason...
        cmd_exec.run_cmd(cmd)
        # remove the .gz suffix from the filename
        reference = reference.with_suffix('')  # Update the reference to remove the .gz suffix

    return reference


def generate_site_tsv(vcf_file: Path, sites_suffix: str, cmd_exec: CommandExecutor = CMD_EXEC) -> Path:
    """This method is a helper to generate a sites file for a given input VCF file

    A simple wrapper around `bcftools query` to generate a .tsv file with columns of CHROM, POS, REF, ALT, AC. This file
    will not contain a header, and thus needs to be considered when reading the file with a :func:`csv.DictReader`.

    :param vcf_file: A Pathlike representation of a VCF or BCF file (can also be gzipped).
    :param sites_suffix: The suffix to append to the sites file.
    :param cmd_exec: An optional alternate CommandExecutor to use for executing system calls. Defaults to the default
        CMD_EXEC from this module.
    :return: A Pathlike to the generated sites file
    """

    sites_file = replace_multi_suffix(vcf_file, sites_suffix)

    cmd = f'bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\n" ' \
          f'-o /test/{sites_file.name} /test/{vcf_file.name}'
    cmd_exec.run_cmd_on_docker(cmd)

    return sites_file


def download_vcf(input_vcf: str, cmd_exec: CommandExecutor = CMD_EXEC) -> Tuple[Path, int]:
    """Helper function that downloads a VCF and it's index, and then calculates the filesize for record-keeping purposes.

    :param input_vcf: An input_vcf file in the form of a DNANexus file ID (file-123456...)
    :param cmd_exec: A CommandExecutor to use for executing system calls on Docker.
    :return: The downloaded VCF file and the size (in bytes) of the VCF file
    """

    # download the VCF files
    vcfpath = InputFileHandler(input_vcf).get_file_handle()

    # if we are using DNA Nexus, find the index file
    if input_vcf == FileType.DNA_NEXUS_FILE:
        find_index(input_vcf, '.tbi')
    # if we are not on DNA Nexus we need to create the index
    else:
        vcfpath_index = vcfpath.with_suffix('.tbi')
        cmd = f'bcftools index -t /test/{vcfpath.name}'
        cmd_exec.run_cmd_on_docker(cmd)

    # Check the file is a vcf.gz
    if not vcfpath.name.endswith('.vcf.gz'):
        raise ValueError(f"File {input_vcf} is not a vcf.gz file!")

    # Check the file is a vcf index
    if not vcfpath_index.name.endswith('.vcf.tbi' or vcfpath_index.name.endswith('.vcf.gz.tbi')):
        raise ValueError(f"File {vcfpath_index.name} is not a tbi file!")

    # Get the size of the VCF file for logging purposes:
    vcf_size = vcfpath.stat().st_size

    return vcfpath, vcf_size


def count_variant_list_and_filter(sites_file: Path, alt_allele_threshold: int) -> Tuple[Path, List[Dict], int, int]:
    """Counts the total number of variants in a given VCF/BCF and filters out sites with excessive alternate alleles

    This method 1st uses the :func:`generate_site_tsv` helper function to query the bcf and print a .tsv file.

    This method then filters sites with excessive alternate alleles. This is because to load them into memory would
    result in excessive memory use and cause the machine running it to crash. This filtered set of variants is then
    written to a sites file to be used for filtering of sites. Filtered sites are retained as a dictionary and
    returned, so they can be written to a file that documents filtered sites. Alleles with a minor allele frequency >
    0.1% are specifically captured in order to know blind spots during subsequent analysis.

    :param sites_file: An *unfiltered* sites file from a VCF
    :param alt_allele_threshold: Number of alternate alleles to allow in a variant before it is excluded
    :param cmd_exec: An optional alternate CommandExecutor to use for executing system calls. Defaults to the default
        CMD_EXEC from this module.
    :return: A tuple of the generated sites file as a Pathlike, a dict of failed sites, the number of variant rows,
        and the number of alternate alleles minus star alleles
    """

    filtered_sites_file = replace_multi_suffix(sites_file, '.sites.txt')

    with sites_file.open('r') as sites_reader, \
            filtered_sites_file.open('w') as filtered_sites_writer:

        filtered_sites_csv = csv.DictWriter(filtered_sites_writer, delimiter='\t', extrasaction='ignore',
                                            fieldnames=['CHROM', 'POS', 'REF', 'ALT'])
        n_vcf_lines = 0
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
                n_var_star = alt_alleles.count('*')

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
                    n_vcf_alternates += (n_var_alt - n_var_star)

        return filtered_sites_file, failed_sites, n_vcf_lines, n_vcf_alternates


def normalise_and_left_correct(vcf_file: Path, site_list: Path, reference_fasta: Path, cmd_exec: CommandExecutor = CMD_EXEC) -> Path:
    """A wrapper for BCFtools norm to left-normalise and split all variants

    Generate a normalised bcf file for all downstream processing using `bcftools norm`:
    -m : splits all multiallelics into separate records
    -f : provides a reference file so bcftools can left-normalise and check records against the reference genome
    -T : restricts the normalisation to a list of sites given by :param site_list:
    --old-rec-tag : sets a tag in the resulting bcf that contains the original record – used for IDing multi-allelics
        after splitting

    :param vcf_file: A list of bcffiles to normalise and left-correct
    :param site_list: A list of sites to restrict to and then normalise and left-correct
    :param reference_fasta: A reference fasta file to use for normalisation and left-correction, output of InputFileHandler
    :param cmd_exec: An optional alternate CommandExecutor to use for executing system calls. Defaults to the default
        CMD_EXEC from this module.
    :return: A Path object to the normalised and left-corrected BCF file
    """

    out_bcf = replace_multi_suffix(vcf_file, '.norm.bcf')
    cmd = f'bcftools norm --threads 8 -w 100 -Ob -m - -f /test/{reference_fasta.name} ' \
          f'-T /test/{site_list.name} ' \
          f'--old-rec-tag MA ' \
          f'-o /test/{out_bcf.name} /test/{vcf_file.name}'
    cmd_exec.run_cmd_on_docker(cmd)

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
    append_last = n_lines % chunk_size < (chunk_size / 2)

    # Then do second iteration to write the actual index files:
    with filtered_sites_file.open('r') as sites_reader:
        # These variables are to control iteration parameters
        current_index_num = 1
        current_end = 0
        write_next_variant = True

        # Store the actual names of each chunk, so we can do the actual bcftools extraction later
        file_chunk_paths = []

        # Collection of files that store different information for i/o:
        # 'sites' is the reader of ALL variants in the BCF
        # 'current_index' stores the actual variants to be extracted for each chunk
        sites = csv.DictReader(sites_reader,
                               delimiter='\t',
                               fieldnames=['chrom', 'pos', 'ref', 'alt', 'ac'],
                               quoting=csv.QUOTE_NONE)

        # Do the iteration itself
        chunk_variant_count = 0
        overall_variant_count = 1
        for site in sites:
            if overall_variant_count == 1:
                current_index_path = replace_multi_suffix(filtered_sites_file,
                                                          f'.chunk{current_index_num}')
                current_index_writer = current_index_path.open('w')
            elif (overall_variant_count - 1) % chunk_size == 0:
                if (n_lines - (overall_variant_count - 1)) < (chunk_size / 2) and append_last is True:
                    logging.info("Appending remaining lines to end of last chunk...")
                else:
                    file_chunk_paths.append(current_index_path)
                    current_index_writer.close()
                    chunk_variant_count = 0
                    current_index_num += 1
                    current_index_path = replace_multi_suffix(filtered_sites_file,
                                                              f'.chunk{current_index_num}')
                    current_index_writer = current_index_path.open('w')

            # Don't write duplicate positions as BCFTools can only process on position rather than position / ref / alt
            if int(site['pos']) == current_end:
                write_next_variant = False

            # Get a running total of all variants EXCEPT star alleles
            if site['alt'] != '*':
                overall_variant_count += 1

            # Allows me to control writing of variants
            if write_next_variant:
                chunk_variant_count += 1
                current_index_writer.write(f"{site['chrom']}\t{site['pos']}\n")
                current_end = int(site['pos'])
            else:
                write_next_variant = True

        # Make sure we close the final writer – in very rare cases the file might be empty due to the last variant
        # being a duplicate position
        current_index_writer.close()
        if chunk_variant_count == 0:
            current_index_path.unlink()
        else:
            file_chunk_paths.append(current_index_path)

    return file_chunk_paths


def split_bcfs(vcf_file: Path, file_chunk_names: List[Path], cmd_exec: CommandExecutor = CMD_EXEC) -> List[Path]:
    """A wrapper around bcftools view that generates smaller chunks from the original VCF file

    We also remove star '*' alleles at this point as we do not use them for any subsequent association testing.

    :param vcf_file: Path to the downloaded VCF to split into chunks
    :param file_chunk_names: A List of individual file chunk Paths of split variant lists
    :param cmd_exec: An optional alternate CommandExecutor to use for executing system calls. Defaults to the default
        CMD_EXEC from this module.
    :return: A List of individual file chunk Pathlikes of split BCF files
    """

    bcf_files = []
    for file_chunk in file_chunk_names:
        out_bcf = file_chunk.with_suffix(f'{file_chunk.suffix}.bcf')
        cmd = f'bcftools view --threads 8 -e "alt==\'*\'" -T /test/{file_chunk.name} ' \
              f'-Ob -o /test/{out_bcf.name} ' \
              f'/test/{vcf_file.name}'
        cmd_exec.run_cmd_on_docker(cmd)
        bcf_files.append(out_bcf)

    return bcf_files


def process_vcf(input_vcf: str, chunk_size: int, alt_allele_threshold: int,
                reference_fasta: Path) -> Tuple[List[dxpy.DXFile], Dict, List[Dict]]:
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
    :param reference_fasta: path to a reference fasta file, output of InputFileHandler
    """

    # 1. Download the VCF
    vcf_path, vcf_size = download_vcf(input_vcf)

    # 2. Generate a list of all variants in the file
    sites_file = generate_site_tsv(vcf_path, '.sites.unfiltered.txt')

    # 3. Filter the site list for large alt size
    orig_sites, failed_sites, n_orig_lines, n_norm_filtered_lines = count_variant_list_and_filter(sites_file,
                                                                                                  alt_allele_threshold)

    # Collate information about this file
    log_info = {'vcf': vcf_path.name, 'dxid': input_vcf, 'n_sites': n_orig_lines,
                'n_final_sites': n_norm_filtered_lines, 'vcf_size': vcf_size}

    # Some WGS-based vcfs are meant to have 0 sites, and we want to capture that here, so we have a full accounting
    if n_norm_filtered_lines == 0:
        final_files = []  # Empty since no splitting will happen
        vcf_path.unlink()

    else:

        # 4. Normalise and left-correct all variants
        norm_bcf = normalise_and_left_correct(vcf_path, orig_sites, reference_fasta)
        vcf_path.unlink()

        norm_sites = generate_site_tsv(norm_bcf, '.norm.sites.txt')

        # 5. Generate reasonable sized (param: chunk_size) lists of variants:
        file_chunk_paths = split_sites(norm_sites, n_norm_filtered_lines, chunk_size)

        # 6. And actually split the files into chunks:
        split_files = split_bcfs(norm_bcf, file_chunk_paths)
        final_files = [generate_linked_dx_file(file) for file in split_files]

    return final_files, log_info, failed_sites


def write_information_files(output_name: Optional[str], n_vcfs: int, infos: List[Dict],
                            skipped_sites: List[List[Dict]], output_dir: Path = Path('./')) -> Tuple[Path, Path]:

    output_name = f'.{output_name}.' if output_name else '.'
    split_info_path = output_dir / f'vcf_info{output_name}tsv'
    skipped_sites_path = output_dir / f'skipped_sites{output_name}tsv'
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

        if len(infos) != len(skipped_sites):
            raise ValueError('Length of infos and skipped_sites do not match!')

        # To future devs: zip does not work here, and I'm too much of a n00b to understand why
        for n_vcf, info in enumerate(infos):
            if info['n_sites'] == 0:
                size_zero_bcf_count += 1
            split_info_csv.writerow(info)
            for site in skipped_sites[n_vcf]:
                site['vcf'] = info['vcf']
            skipped_sites_csv.writerows(skipped_sites[n_vcf])

    LOGGER.info(f'Number of VCFs with 0 sites: {size_zero_bcf_count} '
                f'({(size_zero_bcf_count / n_vcfs) * 100:0.2f}%)')

    return split_info_path, skipped_sites_path


@dxpy.entry_point('main')
def main(input_vcfs: dict, chunk_size: int, alt_allele_threshold: int, output_name: str, human_reference: dict,
         human_reference_index: dict) -> dict:
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
    :return: A dictionary of outputs as specified by dxapp.json
    """

    # Run through each VCF file provided and perform filtering.
    # input_vcfs is a file list of DNANexus file hashes that I dereference below
    # input_vcfs = dxpy.DXFile(input_vcfs['$dnanexus_link'])
    # dxpy.download_dxfile(input_vcfs.get_id(), 'vcf_list.txt')  # Actually download the file

    # download the file
    input_vcfs = InputFileHandler(input_vcfs).get_file_handle()

    reference = ingest_human_reference(human_reference, human_reference_index)

    # Use thread utility to multi-thread this process:
    thread_utility = ThreadUtility(thread_factor=8,
                                   error_message='A splitting thread failed',
                                   incrementor=5)

    # And launch individual jobs
    n_vcfs = 0
    with input_vcfs.open('r') as input_vcf_reader:
        # read in each line of the input file
        for line in input_vcf_reader:
            # Skip the header line if it doesn't start with 'file-' (e.g. the header should be 'vcf / vcf_idx'
            if not line.startswith('file-'):
                continue

            n_vcfs += 1
            # Get the input VCF and index
            input_vcf = line.rstrip().split()[0]

            thread_utility.launch_job(process_vcf,
                                      input_vcf=input_vcf,
                                      chunk_size=chunk_size,
                                      alt_allele_threshold=alt_allele_threshold,
                                      reference_fasta=reference)

    bcf_files = []
    infos = []
    skipped_sites = []
    for result in thread_utility:
        files, info, skipped = result
        bcf_files.extend(files)
        infos.append(info)
        skipped_sites.append(skipped)

    split_info_path, skipped_sites_path = write_information_files(output_name, n_vcfs, infos, skipped_sites)

    # Set output
    output = {'output_vcfs': [dxpy.dxlink(item) for item in bcf_files],
              'run_info': dxpy.dxlink(generate_linked_dx_file(split_info_path)),
              'skipped_sites': dxpy.dxlink(generate_linked_dx_file(skipped_sites_path))}

    return output


if __name__ == '__main__':
    dxpy.run()