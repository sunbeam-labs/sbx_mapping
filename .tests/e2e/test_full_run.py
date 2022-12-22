import csv
import os
import pytest
import subprocess as sp
import tempfile


@pytest.fixture
def dir(pytestconfig):
    return pytestconfig.getoption("dir")


@pytest.fixture
def setup(dir):
    temp_dir = dir if dir else tempfile.mkdtemp()

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


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

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

    output_fp = os.path.join(project_dir, "sunbeam_output/")

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


def test_benchmarks(run_sunbeam):
    (
        human_genome_fp,
        human_copy_genome_fp,
        phix174_genome_fp,
        benchmarks_fp,
    ) = run_sunbeam

    filename = os.listdir(benchmarks_fp)[0]
    with open(os.path.join(benchmarks_fp, filename)) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            assert (
                float(r["cpu_time"]) < 0.5
            ), f"cpu_time for {r['rule']} is higher than 0.5: {r['cpu_time']}"
