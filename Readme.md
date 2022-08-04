# BCFSplitter (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://documentation.dnanexus.com/.

### Table of Contents

- [Introduction](#introduction)
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

### Background

UKBiobank provides a total of 977 individual .vcf.gz files (located at `Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/*.vcf.gz`).
This is the same number as for the 200k release despite the sample size increased by ~2.5x. As such, pipelines designed to
work on the 200k data no longer work as efficiently; we wrote this tool to split provided vcf.gz files into
smaller chunks.

### Dependencies

#### Docker

This applet uses [Docker](https://www.docker.com/) to supply dependencies to the underlying AWS instance
launched by DNANexus. The Dockerfile used to build dependencies is available as part of the MRCEpid organisation at:

https://github.com/mrcepid-rap/dockerimages/blob/main/burdentesting.Dockerfile

This Docker image is built off of the primary 20.04 Ubuntu distribution available via [dockerhub](https://hub.docker.com/layers/ubuntu/library/ubuntu/20.04/images/sha256-644e9b64bee38964c4d39b8f9f241b894c00d71a932b5a20e1e8ee8e06ca0fbd?context=explore).
This image is very light-weight and only provides basic OS installation. Other basic software (e.g. wget, make, and gcc) need
to be installed manually. For more details on how to build a Docker image for use on the UKBiobank RAP, please see:

https://github.com/mrcepid-rap#docker-images

In brief, the primary **bioinformatics software** dependencies required by this Applet (and provided in the associated Docker image)
are:

* [htslib and samtools](http://www.htslib.org/)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)

This list is not exhaustive and does not include dependencies of dependencies and software needed
to acquire other resources (e.g. wget). See the referenced Dockerfile for more information.

#### Resource Files

This applet does not have any external dependencies.

## Methodology

This applet is step 1 (mrc-splitbcf) of the rare variant testing pipeline developed by Eugene Gardner for the UKBiobank
RAP:

![](https://github.com/mrcepid-rap/.github/blob/main/images/RAPPipeline.png)

This applet has three major steps:

1. Generate a list of all variants by chromosome/coordinate in a single vcf file
2. Split this list into smaller chunks (~5000 variants each)
3. Extract variants for each of these chunks into a separate file.

In brief, each vcf file is split into roughly 5000-line chunks. However, a split bcf file can have no fewer than 2500 
sites and no more than 7500 sites. In practice, this results in most original vcf.gz files being split into between 4-5 
smaller .bcfs. Note that this is the _ONLY_ part of the pipeline that I have not parallelized. This was because I did 
this analysis prior to deciding to parallelize all other parts of the codebase.

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

| input      | description                                                       |
|------------|-------------------------------------------------------------------|
| input_vcfs | List of raw input vcf.gz file(s) and their indices from DNA Nexus |

The format of the input_vcfs file is as follows:

```text
file-1234567890ABCDE    file-ABCDE1234567890
```

Where `file-1234567890ABCDE` is the DNANexus file ID for a unprocessed vcf file and `file-ABCDE1234567890` is the 
DNANexus file ID for the unprocessed vcf file's corresponding index.

### Outputs

| output      | description                                                                  |
|-------------|------------------------------------------------------------------------------|
| output_vcfs | All .bcf chunks from the resulting split of all files listed in `input_vcfs` |

output_vcfs is named based on the original file name of the vcf in `input_vcfs` with an additional 'chunk' identifier 
like:

`ukb23157_c1_b0_v1_chunkN.bcf`

### Command line example

If this is your first time running this applet within a project other than "MRC - Variant Filtering", please see our
organisational documentation on how to download and build this app on the DNANexus Research Access Platform:

https://github.com/mrcepid-rap

Running this command is fairly straightforward using the DNANexus SDK toolkit. First generate a list of all VCFs to process
using commands like the following:

```shell
# 1. get a list of the initial VCF
dx ls -l 'Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157*.vcf.gz' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; $F[5] =~ s/.vcf.gz//; print "$F[5]\t$id\n";}' > bcf_list.txt 

# 2. get a concurrent list of the indicies:
dx ls -l 'Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157*.vcf.gz.tbi' | perl -ane 'chomp $_; if ($F[6] =~ /^\((\S+)\)$/) {$id = $1; $F[5] =~ s/.vcf.gz.tbi//; print "$F[5]\t$id\n";}' > idx_list.txt```

# 3. Match the two lists using a custom script (not included in this repo):
matcher.pl -file1 bcf_list.txt -file2 idx_list.txt -r | perl -ane 'chomp $_; print "$F[3]\t$F[1]\n";' > bcfsplitter_input.lst

# 4. Split to a manageble per-size job:
split -l 15 -d ../bcfsplitter_input.lst splitjob_

# 5. Upload to DNAnexus
dx upload splitjob_* --destination batch_files/
```

For the input vcf (provided with the flag `-iinput_vcf`) one can use a file hash from one of the above files:

```shell
dx run mrcepid-bcfsplitter --priority low --destination filtered_vcfs/ \ 
        -iinput_vcfs=file-1234567890ABCDE
```

Brief I/O information can also be retrieved on the command line:

```shell
dx run mrcepid-bcfsplitter --help
```

I have set a sensible (and tested) default for compute resources on DNANexus that is baked into the json used for building 
the app (at `dxapp.json`) so setting an instance type is unnecessary. This current default is for a mem2_ssd1_v2_x64 instance
(64 CPUs, 256 Gb RAM). If necessary to adjust compute resources, one can provide a flag like `--instance-type mem1_ssd1_v2_x4`.

#### Batch Running

This tool is not compatible with batch running at this time.
