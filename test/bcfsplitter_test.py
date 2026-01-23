#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# Prior to using this script PLEASE follow the instructions in the developer readme (Readme.developer.md) carefully.
# This Readme provides instructions on how to regenerate testing data necessary to run these tests.
import csv
import pytest
import shutil
from pathlib import Path

from time import sleep
from typing import Optional, List, Tuple
from pysam import VariantFile

from general_utilities.association_resources import replace_multi_suffix
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor, \
    build_default_command_executor
from bcfsplitter.bcfsplitter import generate_site_tsv, count_variant_list_and_filter, normalise_and_left_correct, \
    split_sites, split_bcfs, write_information_files, ingest_human_reference, download_vcf

test_data_dir = Path(__file__).parent / 'test_data'
CMD_EXEC = build_default_command_executor()

EXPECTED_VCF_VALUES = [{'final': 845, 'missing': 6, 'original': 835,
                        'vcf': test_data_dir / 'test_input1.vcf.gz',
                        'index': test_data_dir / 'test_input1.vcf.gz.tbi'},
                       {'final': 417, 'missing': 32, 'original': 410,
                        'vcf': test_data_dir / 'test_input2.vcf.gz',
                        'index': test_data_dir / 'test_input2.vcf.gz.tbi'}]

test_ref = test_data_dir / 'reference.fasta'
test_fai = test_data_dir / 'reference.fasta.fai'

for vcf_data in EXPECTED_VCF_VALUES:

    assert vcf_data['vcf'].exists()
    assert vcf_data['index'].exists()

assert test_ref.exists()
assert test_fai.exists()

KEEP_TEMP = False  # Set to True if you want to keep temp outputs for debugging


@pytest.fixture(scope='function')
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.
    """
    test_data_source = Path(__file__).parent / "test_data"
    destination = tmp_path / "test_data"
    shutil.copytree(test_data_source, destination)
    monkeypatch.chdir(tmp_path)
    yield tmp_path
    if KEEP_TEMP:
        persistent_dir = Path(__file__).parent / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


def test_ingest_human_reference(temporary_path, cmd_exec=CMD_EXEC):
    """Test the ingestion of a human reference genome.

    :param reference: The name of the reference genome file.
    :param reference_index: The name of the index file for the reference genome.
    """

    reference = temporary_path / 'test_data' / 'reference.fasta'
    reference_index = temporary_path / 'test_data' / 'reference.fasta.fai'

    reference = ingest_human_reference(reference, reference_index, cmd_exec)

    assert reference.exists()
    assert reference_index.exists()


@pytest.mark.parametrize(
    argnames=['input_vcf'],
    argvalues=[
        (test_data_dir / 'test_input1.vcf.gz',),
        (test_data_dir / 'test_input2.vcf.gz',),
    ]
)
def test_download_vcf(temporary_path, input_vcf):
    """Test the download of a VCF file.

    :param input_vcf: The name of the VCF file to be downloaded.
    """

    # Convert to Path objects relative to tmp_data_dir
    input_vcf = temporary_path / input_vcf

    vcfpath, vcf_zie = download_vcf(input_vcf, CMD_EXEC)

    assert vcfpath.exists()
    assert vcf_zie == input_vcf.stat().st_size


@pytest.mark.parametrize(argnames=['sites_suffix', 'vcf_info'],
                         argvalues=zip(['.sites.tsv', '.sites.tsv'], EXPECTED_VCF_VALUES)
                         )
def test_generate_site_tsv(temporary_path, sites_suffix: str, vcf_info: dict, cmd_exec=CMD_EXEC):
    """Test generate_site_tsv

    generate_site_tsv will take a VCF file and generate a site TSV file. This is a simple test to ensure that the
    site TSV file is generated correctly.

    :param tmp_data_dir: A tmp data directory to store the VCF file in.
    :param sites_suffix: The suffix to append to the site TSV file.
    :param vcf_info: A dictionary containing information about the VCF file being tested.
    """

    tmp_vcf = temporary_path / 'test_input.vcf.gz'
    tmp_idx = temporary_path / 'test_input.vcf.gz.tbi'
    vcf_info['vcf'].link_to(tmp_vcf)
    vcf_info['index'].link_to(tmp_idx)

    # Have to do it this way because of temp directories...
    file_list = generate_site_tsv(Path(tmp_vcf), sites_suffix, cmd_exec)
    file_list = tmp_vcf.parent / file_list

    assert file_list.exists()

    assert file_list == replace_multi_suffix(tmp_vcf, '.sites.tsv')
    with file_list.open('r') as test_list:
        for site_num, site in enumerate(test_list):
            data = site.split('\t')
            assert len(data) == 5
            assert int(data[1])

        # Enumerate is 0-based...
        assert site_num + 1 == vcf_info['original']


@pytest.mark.parametrize(argnames=['alt_allele_threshold', 'expected_missing', 'expected_alts'],
                         argvalues=zip([0, 1, 2, 3, 4, 5, 6],
                                       [835, 54, 6, 4, 4, 4, 0],
                                       [0, 781, 845, 851, 851, 851, 875]))
@pytest.mark.parametrize('vcf_info', [EXPECTED_VCF_VALUES[0]])
def test_count_variant_list_and_filter(temporary_path, alt_allele_threshold: int, expected_missing: int,
                                       expected_alts: int, vcf_info: dict, cmd_exec=CMD_EXEC):
    """Test the count_variant_list_and_filter function.

    Here we are primarily testing if the method filters according to alternate allele count properly and returns the
    correct number of sites AFTER filtering.

    Note that we only test one VCF – unnecessary to test both

    :param tmp_data_dir: A tmp data directory to store the VCF file in.
    :param alt_allele_threshold: Max number of alternate alleles to filter on.
    :param expected_missing: Number of sites we expect to be excluded after filtering.
    :param expected_alts: Total number of alternate alleles we expect to be left after filtering. This will be the
        total number of redundant positions in the VCF file.
    :param vcf_info: A dictionary containing information about the VCF file being tested.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    tmp_vcf = temporary_path / 'test_input.vcf.gz'
    tmp_idx = temporary_path / 'test_input.vcf.gz.tbi'
    vcf_info['vcf'].link_to(tmp_vcf)
    vcf_info['index'].link_to(tmp_idx)

    file_list = generate_site_tsv(Path(tmp_vcf), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf.parent / file_list

    filtered_file, failed_sites, n_vcf_lines, \
        n_vcf_alternates = count_variant_list_and_filter(file_list, alt_allele_threshold)
    filtered_file = tmp_vcf.parent / filtered_file
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    assert filtered_file.exists()

    # Check that the number of lines in the filtered file is as expected
    assert n_vcf_lines == vcf_info['original']
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
@pytest.mark.parametrize('vcf_info', [EXPECTED_VCF_VALUES[0]])
def test_normalise_and_left_correct(temporary_path, alt_allele_threshold: int, expected_sites: int,
                                    expected_multi_allelics: int, vcf_info: dict, cmd_exec=CMD_EXEC):
    """Test the normalisation function of BCFtools in the context of this applet.

    :param tmp_data_dir: A tmp data directory to store the VCF file in.
    :param alt_allele_threshold: Max number of alternate alleles to filter on.
    :param expected_sites: Total number of sites AFTER normalisation.
    :param expected_multi_allelics: Expected count of sites with an MA tag after normalisation.
    :param vcf_info: A dictionary containing information about the VCF file being tested.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    tmp_vcf = temporary_path / 'test_input.vcf.gz'
    tmp_idx = temporary_path / 'test_input.vcf.gz.tbi'
    vcf_info['vcf'].link_to(tmp_vcf)
    vcf_info['index'].link_to(tmp_idx)

    # Need to re-run to get the filtered file...
    file_list = generate_site_tsv(Path(tmp_vcf), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf.parent / file_list
    filtered_file, _, _, _ = count_variant_list_and_filter(file_list, alt_allele_threshold)
    normalised_bcf = normalise_and_left_correct(Path(tmp_vcf), Path(filtered_file), test_ref, cmd_exec)
    normalised_bcf = tmp_vcf.parent / normalised_bcf
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    assert normalised_bcf.exists()

    bcf_reader = VariantFile(str(normalised_bcf))
    n_multi_alleleic = 0
    n_ma_fields = 0
    n_var = -1  # Ensure n_var is always defined
    for n_var, rec in enumerate(bcf_reader.fetch()):
        if len(rec.alts) > 1:
            n_multi_alleleic += 1
        if 'MA' in rec.info.keys():
            n_ma_fields += 1

    # n_var is 0-based
    if n_var == -1:
        # No records found
        assert expected_sites == 0
    else:
        assert (n_var + 1) == expected_sites
    assert n_multi_alleleic == 0
    assert n_ma_fields == expected_multi_allelics


@pytest.mark.parametrize(argnames=['chunk_size'],
                         argvalues=zip([50, 100, 500])
                         )
@pytest.mark.parametrize('vcf_info', [EXPECTED_VCF_VALUES[0]])
def test_split_sites_and_split_bcfs(temporary_path, chunk_size: int, vcf_info: dict, cmd_exec=CMD_EXEC):
    """Test splitting of sites and bcf into smaller chunks.

    :param tmp_data_dir: A tmp data directory to store the VCF file in.
    :param chunk_size: The size of the chunk to create
    :param vcf_info: A dictionary containing information about the VCF file being tested.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    tmp_vcf = temporary_path / 'test_input.vcf.gz'
    tmp_idx = temporary_path / 'test_input.vcf.gz.tbi'
    vcf_info['vcf'].link_to(tmp_vcf)
    vcf_info['index'].link_to(tmp_idx)

    file_list = generate_site_tsv(Path(tmp_vcf), '.sites.tsv', cmd_exec)
    file_list = tmp_vcf.parent / file_list
    filtered_file, _, n_lines, norm_lines = count_variant_list_and_filter(file_list, 2)
    normalised_bcf = normalise_and_left_correct(Path(tmp_vcf), Path(filtered_file), test_ref, cmd_exec)
    norm_sites = generate_site_tsv(Path(normalised_bcf), '.norm.sites.txt', cmd_exec)
    norm_sites = tmp_vcf.parent / norm_sites
    chunk_paths = [Path(chunk) for chunk in split_sites(norm_sites, norm_lines, chunk_size)]
    split_files = split_bcfs(Path(normalised_bcf), chunk_paths, cmd_exec)
    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

    expected_sites = vcf_info['final']

    chunk_mod = 0 if expected_sites % chunk_size < (chunk_size / 2) else 1
    expected_chunk_len = (expected_sites // chunk_size) + chunk_mod

    assert len(chunk_paths) == expected_chunk_len
    assert len(split_files) == expected_chunk_len

    num_sites = 0
    for bcf_n, split_bcf in enumerate(split_files):
        chunk_sites = 0
        split_bcf = tmp_vcf.parent / split_bcf
        split_reader = VariantFile(str(split_bcf))
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


@pytest.mark.parametrize(argnames=['output_name', 'vcf_infos'],
                         argvalues=zip(['test_output', None], [EXPECTED_VCF_VALUES] * 2))
def test_write_information_files(temporary_path, output_name: Optional[str], vcf_infos: List[dict], cmd_exec=CMD_EXEC):
    """Test the writing of information files at the conclusion of a run.

    :param tmp_data_dir: A tmp data directory to store the VCF file(s) in.
    :param output_name: The name of the output file. Can be none which will result in a default name being used.
    """

    # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!
    infos = []
    skipped_sites = []
    validation_data = {}
    for vcf_info in vcf_infos:
        tmp_vcf = temporary_path / 'test_input.vcf.gz'
        tmp_idx = temporary_path / 'test_input.vcf.gz.tbi'
        vcf_info['vcf'].link_to(tmp_vcf)
        vcf_info['index'].link_to(tmp_idx)

        file_list = generate_site_tsv(Path(tmp_vcf), '.sites.tsv', cmd_exec)
        file_list = tmp_vcf.parent / file_list

        filtered_file, failed_sites, n_vcf_lines, \
            n_vcf_alternates = count_variant_list_and_filter(file_list, 2)
        # !!! DO NOT MODIFY – THIS IS COMPLEX DUE TO TMP_DIR PATHS !!!

        infos.append({'vcf': vcf_info['vcf'], 'dxid': 'file-1234567890ABCDEFG', 'n_sites': n_vcf_lines,
                     'n_final_sites': n_vcf_alternates, 'vcf_size': tmp_vcf.stat().st_size})
        skipped_sites.append(failed_sites)

        validation_data[vcf_info['vcf']] = {'original': n_vcf_lines, 'final': n_vcf_alternates}

        tmp_vcf.unlink()
        tmp_idx.unlink()
        sleep(2)  # Sleep to ensure the files are deleted before the next test

    info_path, failed_path = write_information_files(output_name, 1, infos, skipped_sites,
                                                     output_dir=temporary_path)

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
        valid_names = [str(x['vcf']) for x in EXPECTED_VCF_VALUES]  # Convert Path to str
        for row in info_csv:
            assert row['vcf'] in valid_names
            assert row['dxid'] == 'file-1234567890ABCDEFG'
            assert int(row['n_sites']) == validation_data[Path(row['vcf'])]['original']  # Use Path for dict lookup
            assert int(row['n_final_sites']) == validation_data[Path(row['vcf'])]['final']

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

        expected_missing = sum([x['missing'] for x in EXPECTED_VCF_VALUES])
        assert total_skipped == expected_missing
