#!/usr/bin/env python
# mrcepid-collapsevariants 0.0.1
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import csv
import dxpy
import math
import subprocess

from typing import List
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor


# This function runs a command on an instance, either with or without calling the docker instance we downloaded
# By default, commands are not run via Docker, but can be changed by setting is_docker = True
from dxpy import DXSearchError
def run_cmd(cmd: str, is_docker: bool = False, stdout_file: str = None, print_cmd = False) -> None:

    # -v here mounts a local directory on an instance (in this case the home dir) to a directory internal to the
    # Docker instance named /test/. This allows us to run commands on files stored on the AWS instance within Docker.
    # This looks slightly different from other versions of this command I have written as I needed to write a custom
    # R script to run STAAR. That means we have multiple mounts here to enable this code to find the script.
    if is_docker:
        cmd = "docker run " \
              "-v /home/dnanexus:/test " \
              "-v /usr/bin/:/prog " \
              "egardner413/mrcepid-burdentesting " + cmd

    if print_cmd:
        print(cmd)

    # Standard python calling external commands protocol
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if stdout_file is not None:
        with open(stdout_file, 'w') as stdout_writer:
            stdout_writer.write(stdout.decode('utf-8'))
        stdout_writer.close()

    # If the command doesn't work, print the error stream and close the AWS instance out with 'dxpy.AppError'
    if proc.returncode != 0:
        print("The following cmd failed:")
        print(cmd)
        print("STDOUT follows\n")
        print(stdout.decode('utf-8'))
        print("STDERR follows\n")
        print(stderr.decode('utf-8'))
        raise dxpy.AppError("Failed to run properly...")


# Utility function to delete files no longer needed from the AWS instance to save space
def purge_file(file: str) -> None:

    cmd = "rm " + file
    run_cmd(cmd)


# This is a helper function to upload a local file and then remove it from the instance.
# This is different than other applets I have written since CADD takes up so much space.
# I don't want to have to use a massive instance costing lots of £s!
def generate_linked_dx_file(file: str) -> dxpy.DXFile:

    linked_file = dxpy.upload_local_file(file)
    purge_file(file)
    return linked_file


# Download required resources for this applet
def ingest_resources() -> None:

    # Bring a prepared docker image into our environment so that we can run commands we need:
    # The Dockerfile to build this image is located at resources/Dockerfile
    cmd = "docker pull egardner413/mrcepid-burdentesting:latest"
    run_cmd(cmd)


# Helper function that downloads a VCF and it's index, then returns the filename prefix
def download_vcf(input_vcf: str) -> str:

    # Unsure if this is safe, but I think is the easiest way to avoid hardcoding a given project ID
    project_ID = dxpy.PROJECT_CONTEXT_ID

    # Set a DX file handler for the VCF file & index chunk
    # Both the DXFile class and download_dxfile method require a project_ID as they directly access UKBB bulk data.
    # The method I have used _hopefully_ avoids hardcoding...
    vcf = dxpy.DXFile(input_vcf,project=project_ID)

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

    cmd = "bcftools query -f \"%CHROM\\t%POS\\n\" -o /test/" + vcfprefix + ".bcf_sites.txt /test/" + vcfprefix + ".vcf.gz"
    run_cmd(cmd, True)


# Essentially is a re-code of the *NIX split command with a few special tweaks to:
# 1. Ensure that BCF files smaller than 2,500 variants do not get created
# 2. Make sure that variants do not get duplicated at the beginning/end of VCFs due to poor handling of MNVs
def split_sites(vcfprefix: str) -> List[str]:

    # Do first iteration to get number of lines – this defines how many files we will write
    sites = csv.DictReader(open(vcfprefix + ".bcf_sites.txt", "r"), delimiter="\t", fieldnames = ["chrom","pos"], quoting = csv.QUOTE_NONE)
    n_lines = 0
    for site in sites:
        n_lines += 1
    if n_lines % 5000 < 2500:
        append_last = True
    else:
        append_last = False

    # Then do second iteration to write the actual index files:
    # These variables are to control iteration parameters
    current_index_num = 1
    current_end = 0
    write_next_variant = True

    # Store the actual names of each chunk, so we can do the actual bcftools extraction later
    file_chunk_names = [vcfprefix + "_chunk" + str(current_index_num)]

    # Collection of files that store different information for i/o
    sites = csv.DictReader(open(vcfprefix + ".bcf_sites.txt", "r"),  # 'sites' is the reader of ALL variants in the BCF
                           delimiter="\t",
                           fieldnames=["chrom", "pos"],
                           quoting=csv.QUOTE_NONE)
    current_index = open(vcfprefix + "_chunk" + str(current_index_num), "w")  # Stores the actual variants to be extracted for each chunk

    # Do the iteration itself
    for site in sites:
        if sites.line_num == 1:
            current_start = site['pos']
        elif (sites.line_num - 1) % 5000 == 0:
            if (n_lines - sites.line_num) < 2500 and append_last is True:
                print("Appending remaining lines to end of last chunk...")
            else:
                current_index.close()
                current_index_num += 1
                current_index = open(vcfprefix + "_chunk" + str(current_index_num), "w")
                file_chunk_names.append(vcfprefix + "_chunk" + str(current_index_num))
                if current_end == site['pos']:
                    write_next_variant = False
                current_start = site['pos']
        # Allows me to control writing of variants
        if write_next_variant:
            current_index.write("%s\t%s\n" % (site['chrom'], site['pos']))
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
        cmd = "bcftools view --threads 4 -T /test/" + file + " -Ob -o /test/" + vcfprefix + "_chunk" + str(current_chunk) + ".bcf /test/" + vcfprefix + ".vcf.gz"
        run_cmd(cmd, True)
        bcf_files.append(generate_linked_dx_file(vcfprefix + "_chunk" + str(current_chunk) + ".bcf"))
        current_chunk += 1

    purge_file(vcfprefix + ".vcf.gz")
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

    # Get threads available to this instance
    threads = os.cpu_count()
    print('Number of threads available: %i' % threads)

    # Ingest requisite resources:
    ingest_resources()

    # Run through each VCF file provided and perform filtering.
    # input_vcfs is simple a file list of DNANexus file hashes that I dereference below
    input_vcfs = dxpy.DXFile(input_vcfs)
    dxpy.download_dxfile(input_vcfs.get_id(), "vcf_list.txt")  # Actually download the file
    input_vcf_reader = open("vcf_list.txt", 'r')

    # Now build a thread worker that contains as many threads, divided by 4 that have been requested since each bcftools
    # instance takes 4 threads and 1 thread for monitoring
    available_workers = math.floor((threads - 1) / 4)
    executor = ThreadPoolExecutor(max_workers=available_workers)

    # And launch the requested threads
    future_pool = []
    for line in input_vcf_reader:
        input_vcf = line.rstrip()

        # Perform filtering/annotation in a separate thread for each VCF file
        future_pool.append(executor.submit(process_vcf,
                                           input_vcf=input_vcf))

    input_vcf_reader.close()
    print("All threads submitted...")

    bcf_files = []
    for future in futures.as_completed(future_pool):
        try:
            result = future.result()
            bcf_files.extend(result)
        except Exception as err:
            print("A thread failed...")
            print(Exception, err)

    # Set output
    output = {"output_vcfs": [dxpy.dxlink(item) for item in bcf_files]}

    return output


dxpy.run()
