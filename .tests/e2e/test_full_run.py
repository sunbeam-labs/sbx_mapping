import os
import pytest
import shutil
import subprocess as sp
import sys
import tempfile


@pytest.fixture
def setup():
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

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    output_fp = os.path.join(project_dir, "sunbeam_output")

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "all_mapping",
                "--directory",
                temp_dir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
        sys.exit(e)

    shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
    shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")

    human_genome_fp = os.path.join(output_fp, "mapping/human/")
    human_copy_genome_fp = os.path.join(output_fp, "mapping/human_copy/")
    phix174_genome_fp = os.path.join(output_fp, "mapping/phix174/")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield human_genome_fp, human_copy_genome_fp, phix174_genome_fp, benchmarks_fp


@pytest.fixture
def expected_file_list():
    yield sorted(
        [
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
        ]
    )


def test_full_run(run_sunbeam, expected_file_list):
    (
        human_genome_fp,
        human_copy_genome_fp,
        phix174_genome_fp,
        benchmarks_fp,
    ) = run_sunbeam
    output_files = expected_file_list

    # Check output
    assert sorted(os.listdir(human_genome_fp)) == output_files
    assert sorted(os.listdir(human_copy_genome_fp)) == output_files
    assert sorted(os.listdir(phix174_genome_fp)) == output_files
