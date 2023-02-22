#!/usr/bin/env python
# mrcepid-collapsevariants 0.0.1
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import csv
import logging
import re

import dxpy

from pathlib import Path
from typing import List

# Download required resources for this applet
from general_utilities.association_resources import run_cmd, generate_linked_dx_file
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()


def ingest_resources() -> None:

    # Bring a prepared docker image into our environment so that we can run commands we need:
    # The Dockerfile to build this image is located at resources/Dockerfile
    cmd = "docker pull egardner413/mrcepid-burdentesting:latest"
    run_cmd(cmd, is_docker=False)


# Helper function that downloads a VCF and it's index, then returns the filename prefix
def download_vcf(input_vcf: str) -> str:

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
    dxpy.download_dxfile(dxid=vcf.get_id(), project=project_ID, filename=vcf.describe()['name'])
    dxpy.download_dxfile(dxid=vcfidx.get_id(), project=project_ID, filename=vcfidx.describe()['name'])

    # Set a prefix name for all files so that we can output a standard-named file:
    vcfprefix = vcf.describe()['name'].split(".vcf.gz")[0]

    return vcfprefix


# Just generates a txt file of all variants in the provided BCF file
def generate_variant_list(vcfprefix: str) -> None:

    cmd = f'bcftools query -f "%CHROM\\t%POS\\n" -o /test/{vcfprefix}.bcf_sites.txt /test/{vcfprefix}.vcf.gz'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')


# Essentially is a re-code of the *NIX split command with a few special tweaks to:
# 1. Ensure that BCF files smaller than 2,500 variants do not get created
# 2. Make sure that variants do not get duplicated at the beginning/end of VCFs due to poor handling of MNVs
def split_sites(vcfprefix: str) -> List[str]:

    # Run wc -l to get number of lines – this defines how many files we will write
    run_cmd(f'wc -l {vcfprefix}.bcf_sites.txt', is_docker=False, stdout_file=f'wc_{vcfprefix}.txt')
    with Path(f'wc_{vcfprefix}.txt').open('r') as sites_reader:
        n_lines = -1
        for line in sites_reader:
            line_count_match = re.match('^\s*(\\d+)\s', line)
            if line_count_match:
                n_lines = int(line_count_match.group(1))
        if n_lines == -1:
            raise ValueError(f'could not determine file length of {vcfprefix}.bcf_sites.txt.')

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


# Just a wrapper around bcftools that generates the actual chunks for each VCF file
def split_bcfs(vcfprefix: str, file_chunk_names: List[str]) -> List[dxpy.DXFile]:

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


# Helper function that enables multithreading in this applet
def process_vcf(input_vcf: str) -> List[dxpy.DXFile]:

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
def main(input_vcfs: dict) -> dict:

    # Ingest requisite resources:
    ingest_resources()

    # Run through each VCF file provided and perform filtering.
    # input_vcfs is simple a file list of DNANexus file hashes that I dereference below
    input_vcfs = dxpy.DXFile(input_vcfs)
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


dxpy.run()
