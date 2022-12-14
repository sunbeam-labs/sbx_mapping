import os
import pytest
import shutil
import subprocess as sp
import tempfile


@pytest.fixture
def setUp():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")
    genomes_fp = os.path.abspath(".tests/data/hosts/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

    config_str = f"sbx_mapping: {{genomes_fp: {genomes_fp}}}"

    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    output_fp = os.path.join(project_dir, "sunbeam_output")
    # shutil.copytree(".tests/data/sunbeam_output", output_fp)

    human_genome_fp = os.path.join(output_fp, "mapping/human/")
    human_copy_genome_fp = os.path.join(output_fp, "mapping/human_copy/")
    phix174_genome_fp = os.path.join(output_fp, "mapping/phix174/")

    yield temp_dir, project_dir, human_genome_fp, human_copy_genome_fp, phix174_genome_fp

    shutil.rmtree(temp_dir)


@pytest.fixture
def expected_file_list():
    yield [
        "random.raw.bcf",
        "dummybfragilis.bam.bai",
        "random.bam",
        "dummybfragilis.raw.bcf",
        "dummybfragilis.bam",
        "dummyecoli.raw.bcf",
        "dummyecoli.bam.bai",
        "dummyecoli.bam",
        "random.bam.bai",
        "coverage.csv",
    ].sort()


def test_full_run(setup, expected_file_list):
    (
        temp_dir,
        project_dir,
        human_genome_fp,
        human_copy_genome_fp,
        phix174_genome_fp,
    ) = setup
    output_files = expected_file_list

    # Run the test job
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_mapping",
            "--directory",
            temp_dir,
        ]
    )

    # Check output
    assert os.listdir(human_genome_fp).sort() == output_files
    assert os.listdir(human_copy_genome_fp).sort() == output_files
    assert os.listdir(phix174_genome_fp).sort() == output_files
