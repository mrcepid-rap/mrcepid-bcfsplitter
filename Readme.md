# BCFSplitter (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
  * [Changelog](#changelog)
  * [Background](#background)
  * [Dependencies](#dependencies)
    + [Docker](#docker)
    + [Resource Files](#resource-files)
- [Methodology](#methodology)
- [Running on DNANexus](#running-on-dnanexus)
  * [Inputs](#inputs)
  * [Outputs](#outputs)
  * [Command line example](#command-line-example)
    + [Batch Running](#batch-running)


## Introduction

This applet splits raw vcf.gz files provided by the UKBiobank platform into more manageable chunks for downstream 
processing.

This README makes use of DNANexus file and project naming conventions. Where applicable, an object available on the DNANexus
platform has a hash ID like:

* file – `file-1234567890ABCDEFGHIJKLMN`
* project – `project-1234567890ABCDEFGHIJKLMN`

Information about files and projects can be queried using the `dx describe` tool native to the DNANexus SDK:

```shell
dx describe file-1234567890ABCDEFGHIJKLMN
```

**Note:** This README pertains to data included as part of the DNANexus project "MRC - Variant Filtering" (project-G2XK5zjJXk83yZ598Z7BpGPk)

### Changelog

* v2.0.0
  * Added support for WGS data. This comes with several changes under-the-hood to support issues with the larger WGS data. These changes will not effect processing of WES data, but please see below for specific changes:
    * Modified standard number of CPUs per-vcf to 8. We have done extensive testing to ensure that the process is $O(n^2)$ in terms of CPU usage – e.g., 2 CPUs = 60m / VCF, 4 CPUs = 30m / VCF
    * Fixed an issue with sites with large number of alternate alleles.
      * Issue can lead to excessive memory requirements when processing some VCFs.
      * We have added a method to filter sites with excessive alternate alleles, which can be triggered using the `alt_allele_threshold` parameter.
      * This parameter is set to 99999 by default, so will not affect WES processing
      * We recommend a setting of 15 when analysing WGS data
      * All sites that are filtered will be reported in the `skipped_sites.{output_name}.tsv`
  * Several quality-of-life changes to the codebase to improve readability and maintainability
  * Implemented newer version of `general_utilities`
  * Added docstrings to all functions

* v1.0.0
  * Initial numbered release, see git logs for previous changes.

### Background

UKBiobank provides individual .vcf.gz files. This is the same number as for the 200k release despite the sample size 
increased by ~2.5x. As such, pipelines designed to work on the 200k data no longer work as efficiently; we wrote 
this tool to split provided vcf.gz files into smaller chunks. Default file locations on the RAP are as follows:

* WES: `/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/`
* WGS: `/Bulk/DRAGEN WGS/Whole genome variant call files (GVCFs) (DRAGEN) [500k release]/`

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/burdentesting.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

or 

https://github.com/mrcepid-rap/dockerimages

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [htslib and samtools](http://www.htslib.org/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)

This list is not exhaustive and does not include dependencies of dependencies and software needed to acquire other 
resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

This applet does not have any external dependencies.

## Methodology

This applet is step 1 (mrcepid-splitbcf) of the rare variant testing pipeline developed by Eugene Gardner for the 
UKBiobank RAP:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.v3.png)

This applet has three major steps:

1. Generate a list of all variants by chromosome/coordinate in a single vcf file
   2. If required & requested, filter out sites with excessive alternate alleles
2. Normalise and left correct VCFs, while removing sites with excessive alternate alleles if requested
3. Split this list into smaller chunks, by default (see `chunk_size` parameter) ~5000 variants each
4. Extract variants for each of these chunks into a separate file.

In brief, each vcf file is split into roughly 5000-line chunks by default. However, a split bcf file cannot have 
< `chunk_size` / 2 sites or > `chunk_size` + (`chunk_size` / 2) sites. In practice, and for WES data, this results in
most original vcf.gz files being split into between 4-5 smaller .bcfs.

**BIG Note** There is a fix implemented in this pipeline that solves a previous bug that would duplicate a site if both
of the following were true:

1. The site has an identical position value to another site in the vcf (i.e. split multiallelics)
2. The site and it's position duplicate **SPAN** the junction of a bcf split
    * e.g. if variant 1 of a multiallelic pair is variant no. 5000 and variant 2 is variant no. 5001

This bug resulted in both variants in the pair being present in both split bcf files. In theory, this could also 
happen to a 3 variant multiallelic, but I have not encountered such a situation. As noted, I have since fixed this issue 
in the resulting VCFs using a custom script. I am leaving this documentation in here for posterity and as a warning in 
case this error is ever encountered again.

## Running on DNANexus

### Inputs

| input                 | Optional | Default                                                        | description                                                                                                                                           |
|-----------------------|----------|----------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| input_vcfs            | False    | N/A                                                            | List of raw input vcf.gz file(s) from DNA Nexus                                                                                                       |
| chunk_size            | True     | 5000                                                           | The number of variants to include per-output BCF produced by this applet.                                                                             |
| output_name           | True     | None                                                           | Additional string to use when creating information files (run info, skipped_sites). This does NOT modify the name of split VCF files.                 |
| alt_allele_threshold  | True     | 99,999                                                         | The maximum number of alternate alleles in a single variant before filtering a site. Uses '>' (greater than) to determine threshold.                  |
| human_reference       | True     | project-Fx2x0fQJ06KfqV7Y3fFZq1jp:file-Fx2x270Jx0j17zkb3kbBf6q2 | dxfile / path pointing to the reference genome the provided VCFs are aligned to                                                                       |
| human_reference_index | True     | project-Fx2x0fQJ06KfqV7Y3fFZq1jp:file-Fx2x21QJ06f47gV73kZPjkQQ | dxfile / path pointing to the reference genome index (.fai) the provided VCFs are aligned to                                                          |
| testing_script        | True     | None                                                           | Invoke the bcfsplitter test suite by providing a script compatible with the 'pytest' module. DO NOT use this flag unless you know what you are doing! |
| testing_directory     | True     | None                                                           | Directory name containing test files. DO NOT use this flag unless you know what you are doing!                                                        |

The format of the input_vcfs file is as follows:

```text
file-1234567890ABCDEFGHIJ
file-ABCDEFGHIJ1234567890
file-0987654321JIHGFEDCBA
file-JIHGFEDCBA0987654321
```

Where `file-1234567890ABCDEEFGHIJ` is the DNANexus file ID for an unprocessed vcf file.

#### Generating an input_vcfs List

Please find an example (for WES data) of generating the input list(s) for the `input_vcfs` parameter below:

```shell
# 1. get a list of all VCFs
dx ls -l 'Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157*.vcf.gz' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; $F[5] =~ s/.vcf.gz//; print "$F[5]\t$id\n";}' > bcf_list.txt 

# 2. Split to a manageble per-size job:
split -l 18 -d ../bcfsplitter_input.lst splitjob_

# 3. Upload to DNAnexus
dx upload splitjob_* --destination batch_files/
```

**Note**: `-l 18` is used above to ensure that all CPUs on the default instance type are used. Please see 
[VM Resources](#vm-resources) for more information on why this parameter was used as an example.

#### Alt Allele Threshold

The `alt_allele_threshold` parameter is used to filter out sites with excessive alternate alleles. When using this tool 
with WGS data, we suggest using a value of `15`. We have optimised this parameter to allow for the applet to function 
on a `mem1_ssd1` instance type (`c5d` family on AWS). This is to use low memory machines, which are not typicaly in
[demand on AWS](https://aws.amazon.com/ec2/spot/instance-advisor/). Thus, one can often use low priority on DNANexus
to save costs.

### Outputs

| output        | description                                                                                                                      |
|---------------|----------------------------------------------------------------------------------------------------------------------------------|
| output_vcfs   | All .bcf chunks from the resulting split of all files listed in `input_vcfs`                                                     |
| run_info      | Summary statistics for the each VCF file split as part of this process. Will be named like `vcf_info.{output_name}.tsv`          |
| skipped_sites | A .tsv file containing sites with > alt_allele_threshold alternate alleles. Will be named like `skipped_sites.{output_name}.tsv` |

output_vcfs is named based on the original file name of the vcf in `input_vcfs` with an additional 'chunk' identifier 
like:

`ukb23157_c1_b0_v1_chunkN.bcf`

run_info contains the following columns:

    vcf: The original vcf file
    dxid: The dxid of the original vcf file
    n_sites: Number of sites found in the original file
    n_final_sites: Number of sites after filtering for alt_allele_threshold
    norm_n_alts: Number of alternate alleles found in the normalised file – this SHOULD equal n_norm_sites
    vcf_size: Size of the original vcf file in bytes

skipped_sites contains the following columns:
  
    vcf: The original vcf file
    chrom: Chromosome of the site
    pos: Position of the site
    alts: Alternate alleles with MAF ~> 0.1% as a comma-separated list
    acs: Allele counts for each alternate allele listed in 'alts' as a comma-separated list

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

For the input vcf (provided with the flag `-iinput_vcf`) one can use a file hash from one of the above files:

```shell
# 1. WES example: 
dx run mrcepid-bcfsplitter --priority low --destination filtered_vcfs/ \ 
        -iinput_vcfs=file-1234567890ABCDE -ioutput_name=list1
        
# 2. WGS example:
dx run mrcepid-bcfsplitter --priority low --destination filtered_vcfs/ \ 
        -iinput_vcfs=file-1234567890ABCDE -ioutput_name=list1 -ialt_allele_threshold=15 -ichunk_size=2500
```

Brief I/O information can also be retrieved on the command line:

```shell
dx run mrcepid-bcfsplitter --help
```

### VM Resources

Sensible (and tested) defaults for compute resources on DNANexus that is baked into the json used for building 
the app (at `dxapp.json`). The current default is for a mem1_ssd1_v2_x72 instance (72 CPUs, 144 Gb RAM). This machine
is typically not in high demand on AWS and should allow for jobs to run on low priority. If this machine does 
appear to be in high demand, switch to another instance type in the c5d family of machines (mem1_ssd1 on DNANexus).

Each VCF requires 8 CPUs to run. Thus, when determining the number of VCFs to run per job, one should take the number
of CPUs available on the machine and divide by 8. For example, since the default machine has 72 CPUs, one should run a 
multiple of 9 VCFs per job (e.g., 9, 18, 27, 36, ...). Note that runtime is roughly linear. Thus, if you have 18 VCFs, 
you should expect the job to take double the time of a job with 9 VCFs. This is important to consider when using low 
priority as the longer a job runs, the more likely it is to be interrupted. Using the default instance type, we have 
had good results on low priority with 18 VCFs / job. This strikes a reasonable balance between number of jobs to monitor
(DNANexus is limited to 100 jobs / user) and the likelihood of using low priority.

If necessary to adjust compute resources, one can provide a flag like `--instance-type mem1_ssd1_v2_x36`.
