#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# Prior to using this script PLEASE follow the instructions in the developer readme (Readme.developer.md) carefully.
# This Readme provides instructions on how to regenerate testing data necessary to run these tests.
import csv
import pytest

from pathlib import Path
from typing import Optional
from pysam import VariantFile

from general_utilities.job_management.command_executor import DockerMount, CommandExecutor
from bcfsplitter.bcfsplitter import generate_site_tsv, count_variant_list_and_filter, normalise_and_left_correct, \
    split_sites, split_bcfs, write_information_files

test_data_dir = Path(__file__).parent / 'test_data'

test_vcf = test_data_dir / 'test_input.vcf.gz'
test_idx = test_data_dir / 'test_input.vcf.gz.tbi'
test_ref = test_data_dir / 'reference.fasta'
test_fai = test_data_dir / 'reference.fasta.fai'
test_vcf_len = 835

assert test_vcf.exists()
assert test_idx.exists()
assert test_ref.exists()
assert test_fai.exists()


@pytest.fixture(scope='function')
def tmp_vcf_file(tmp_path_factory) -> Path:
    """Is a fixture to move testing files into a tmp_dir so tests do not clutter the test_data directory.

    This fixture just symlinks the VCF file in `test_data/` to a tmp_dir.

    :param tmp_path_factory: A tmp_path_factory object from pytest
    :return: A Pathlike to the temporary VCF file.
    """

    tmp_data = tmp_path_factory.mktemp('data')
    tmp_vcf = tmp_data / test_vcf.name
    tmp_idx = tmp_data / test_idx.name
    tmp_ref = tmp_data / test_ref.name
    tmp_fai = tmp_data / test_fai.name
    test_vcf.link_to(tmp_vcf)
    test_idx.link_to(tmp_idx)
    test_ref.link_to(tmp_ref)
    test_fai.link_to(tmp_fai)

    return tmp_vcf


@pytest.mark.parametrize(argnames=['sites_suffix'],
                         argvalues=zip(['.sites.tsv'])
                         )
def test_generate_site_tsv(tmp_vcf_file, sites_suffix: str):
    """Test generate_site_tsv

    generate_site_tsv will take a VCF file and generate a site TSV file. This is a simple test to ensure that the
    site TSV file is generated correctly.

    :param tmp_vcf_file: A sym-link to the VCF file to generate the site TSV file from.
    :param sites_suffix: The suffix to append to the site TSV file.
    """

    test_mount = DockerMount(tmp_vcf_file.parent, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # Have to do it this way because of temp directories...
    file_list = generate_site_tsv(Path(tmp_vcf_file.name), sites_suffix, cmd_exec)
    file_list = tmp_vcf_file.parent / file_list

    assert file_list.exists()
    assert file_list.name == 'test_input.sites.tsv'
    with file_list.open('r') as test_list:
        for site_num, site in enumerate(test_list):
            data = site.split('\t')
            assert len(data) == 5
            assert int(data[1])

        # Enumerate is 0-based...
        assert site_num + 1 == test_vcf_len


@pytest.mark.parametrize(argnames=['alt_allele_threshold', 'expected_missing', 'expected_alts'],
                         argvalues=zip([0, 1, 2, 3, 4, 5, 6],
                                       [835, 54, 6, 4, 4, 4, 0],
                                       [0, 781, 845, 851, 851, 851, 875]))
def test_count_variant_list_and_filter(tmp_vcf_file, alt_allele_threshold: int, expected_missing: int,
                                       expected_alts: int):
    """Test the count_variant_list_and_filter function.

    Here we are priarily testing if the method filters according to alternate allele count properly and returns the
    correct number of sites AFTER filtering.

    :param tmp_vcf_file: A sym-link to the VCF file to generate the site TSV file from.
    :param alt_allele_threshold: Max number of alternate alleles to filter on.
    :param expected_missing: Number of sites we expect to be excluded after filtering.
    :param expected_alts: Total number of alternate alleles we expect to be left after filtering. This will be the
        total number of redundant positions in the VCF file.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    test_mount = DockerMount(tmp_vcf_file.parent, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    file_list = generate_site_tsv(Path(tmp_vcf_file.name), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf_file.parent / file_list

    filtered_file, failed_sites, n_vcf_lines, \
        n_vcf_alternates = count_variant_list_and_filter(file_list, alt_allele_threshold)
    filtered_file = tmp_vcf_file.parent / filtered_file
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    assert filtered_file.exists()

    # Check that the number of lines in the filtered file is as expected
    assert n_vcf_lines == test_vcf_len
    assert n_vcf_alternates == expected_alts

    # Check the catalogued failed sites are as expected
    assert len(failed_sites) == expected_missing
    for site in failed_sites:
        assert sorted(set(site.keys())) == sorted({'chrom', 'pos', 'alts', 'acs'})
        assert int(site['pos'])
        alts = site['alts'].split(',')
        counts = site['acs'].split(',')
        assert len(alts) == len(counts)
        if alts[0] != '':
            counts = list(map(int, counts))
            for n, alt in enumerate(alts):
                assert counts[n] > 1000


@pytest.mark.parametrize(argnames=['alt_allele_threshold', 'expected_sites', 'expected_multi_allelics'],
                         argvalues=zip([1, 2, 3, 4, 5, 6], [781, 877, 883, 883, 883, 907], [0, 96, 102, 102, 102, 126]))
def test_normalise_and_left_correct(tmp_vcf_file, alt_allele_threshold: int, expected_sites: int,
                                    expected_multi_allelics: int):
    """Test the normalisation function of BCFtools in the context of this applet.

    :param tmp_vcf_file: A sym-link to the VCF file to generate the site TSV file from.
    :param alt_allele_threshold: Max number of alternate alleles to filter on.
    :param expected_sites: Total number of sites AFTER normalisation.
    :param expected_multi_allelics: Expected count of sites with an MA tag after normalisation.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    test_mount = DockerMount(tmp_vcf_file.parent, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # Need to re-run to get the filtered file...
    file_list = generate_site_tsv(Path(tmp_vcf_file.name), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf_file.parent / file_list
    filtered_file, _, _, _ = count_variant_list_and_filter(file_list, alt_allele_threshold)
    normalised_bcf = normalise_and_left_correct(Path(tmp_vcf_file.name), Path(filtered_file.name), cmd_exec)
    normalised_bcf = tmp_vcf_file.parent / normalised_bcf
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    assert normalised_bcf.exists()

    bcf_reader = VariantFile(normalised_bcf)
    n_multi_alleleic = 0
    n_ma_fields = 0
    for n_var, rec in enumerate(bcf_reader.fetch()):
        if len(rec.alts) > 1:
            n_multi_alleleic += 1
        if 'MA' in rec.info.keys():
            n_ma_fields += 1

    # n_var is 0-based
    assert (n_var + 1) == expected_sites
    assert n_multi_alleleic == 0
    assert n_ma_fields == expected_multi_allelics


@pytest.mark.parametrize(argnames=['chunk_size', 'expected_sites'],
                         argvalues=zip([50, 100, 500],
                                       [845, 845, 845])
                         )
def test_split_sites_and_split_bcfs(tmp_vcf_file, chunk_size: int, expected_sites: int):
    """Test splitting of sites and bcf into smaller chunks.

    :param tmp_vcf_file: A sym-link to the VCF file to generate the site TSV file from.
    :param chunk_size: The size of the chunk to create
    :param expected_sites: Number of expected sites across all split chunks. This is always the same but parameterised
        for future development if required.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    test_mount = DockerMount(tmp_vcf_file.parent, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    file_list = generate_site_tsv(Path(tmp_vcf_file.name), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf_file.parent / file_list
    filtered_file, _, n_lines, norm_lines = count_variant_list_and_filter(file_list, 2)
    normalised_bcf = normalise_and_left_correct(Path(tmp_vcf_file.name), Path(filtered_file.name), cmd_exec)
    norm_sites = generate_site_tsv(Path(normalised_bcf.name), '.norm.sites.txt', cmd_exec)
    norm_sites = tmp_vcf_file.parent / norm_sites
    chunk_paths = [Path(chunk.name) for chunk in split_sites(norm_sites, norm_lines, chunk_size)]
    split_files = split_bcfs(Path(normalised_bcf.name), chunk_paths, cmd_exec)
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    chunk_mod = 0 if expected_sites % chunk_size < (chunk_size / 2) else 1
    expected_chunk_len = (expected_sites // chunk_size) + chunk_mod

    assert len(chunk_paths) == expected_chunk_len
    assert len(split_files) == expected_chunk_len

    num_sites = 0
    for bcf_n, split_bcf in enumerate(split_files):
        chunk_sites = 0
        split_bcf = tmp_vcf_file.parent / split_bcf
        split_reader = VariantFile(split_bcf)

        for rec in split_reader.fetch():
            assert '*' not in rec.alts
            chunk_sites += 1
            num_sites += 1

        if bcf_n == len(split_files) - 1:
            if chunk_mod == 1:
                assert chunk_sites < chunk_size
            else:
                assert chunk_sites > chunk_size
        else:
            assert chunk_sites == chunk_size

    assert num_sites == expected_sites


@pytest.mark.parametrize(argnames=['output_name', 'expected_sites', 'expected_missing'],
                         argvalues=zip(['test_output', None], [845, 845], [6, 6])
                         )
def test_write_information_files(tmp_vcf_file, output_name: Optional[str], expected_sites: int, expected_missing: int):
    """Test the writing of information files at the conclusion of a run.

    :param tmp_vcf_file: A sym-link to the VCF file to generate the site TSV file from.
    :param output_name: The name of the output file. Can be none which will result in a default name being used.
    :param expected_sites: The number of sites we expect to be in the final VCF file.
    :param expected_missing: The number of sites we expect to be missing from the final VCF file but listed in the
        failed sites file.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    test_mount = DockerMount(tmp_vcf_file.parent, Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    file_list = generate_site_tsv(Path(tmp_vcf_file.name), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf_file.parent / file_list

    filtered_file, failed_sites, n_vcf_lines, \
        n_vcf_alternates = count_variant_list_and_filter(file_list, 2)
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    info = [{'vcf': tmp_vcf_file.name, 'dxid': 'file-1234567890ABCDEFG', 'n_sites': n_vcf_lines,
             'n_final_sites': n_vcf_alternates, 'vcf_size': tmp_vcf_file.stat().st_size}]

    info_path, failed_path = write_information_files(output_name, 1, info, failed_sites,
                                                     output_dir=tmp_vcf_file.parent)

    assert info_path.exists()
    assert failed_path.exists()

    if output_name:
        assert info_path.name == f'vcf_info.{output_name}.tsv'
        assert failed_path.name == f'skipped_sites.{output_name}.tsv'
    else:
        assert info_path.name == f'vcf_info.tsv'
        assert failed_path.name == f'skipped_sites.tsv'

    with info_path.open('r') as info_reader:
        info_csv = csv.DictReader(info_reader, delimiter='\t')
        for row in info_csv:
            assert row['vcf'] == tmp_vcf_file.name
            assert row['dxid'] == 'file-1234567890ABCDEFG'
            assert int(row['n_sites']) == test_vcf_len
            assert int(row['n_final_sites']) == expected_sites

    with failed_path.open('r') as failed_reader:
        failed_csv = csv.DictReader(failed_reader, delimiter='\t')
        total_skipped = 0
        for row in failed_csv:
            total_skipped += 1

            alts = row['alts'].split(',')
            counts = row['acs'].split(',')
            assert len(alts) == len(counts)
            if alts[0] != '':
                counts = list(map(int, counts))
                for n, alt in enumerate(alts):
                    assert counts[n] > 1000

        assert total_skipped == expected_missing
