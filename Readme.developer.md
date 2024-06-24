# mrcepid-filterbcf Developer Readme

What this tests: Function-level tests for the `bcfsplitter` applet.

What this does not test: End-to-end tests for the `bcfsplitter` applet. Always make sure to run the applet on the DNAnexus platform before deploying it!

## Test Data Generation

We generated a test VCF file for this applet from the 1000 Genomes Phase 3 Data. We used the following commands:

```bash
bcftools view -Oz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-0
5_chr7.recalibrated_variants.vcf.gz "chr7:100679507-100694250" > test_input.vcf.gz
```

**Note**: To avoid clogging the test directory with temporary files, we use a tmp_dir from pytest. This means that file 
paths in the testing script are a bit _weird_. Do not modify these paths unless you know what you are doing or ask the 
developers.

## Required external files

Tests require a human reference genome. This file is too big to store in the repo. Please download it from DNANexus / external site:

```bash
cd test/test_data/
# .fa
dx download file-Fx2x270Jx0j17zkb3kbBf6q2
# .fai
dx download file-Fx2x21QJ06f47gV73kZPjkQQ

mv hs38DH.fa.gz reference.fasta.gz 
gunzip reference.fasta.gz 
mv hs38DH.fa.fai reference.fasta.fai

```

## Running tests

I use pycharm to run tests through a GUI – but if required to run on the command line, make sure you are in the `test/' 
directory and run:

```bash
pytest bcfsplitter_test.py
```

You will need to ensure the `bcfsplitter` module and `general_utilities` are properly installed and in your python path. 