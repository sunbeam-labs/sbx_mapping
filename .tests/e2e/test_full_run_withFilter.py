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
    hosts_fp = os.path.abspath(".tests/data/hosts/")
    genomes_fp = os.path.abspath(".tests/data/ref/")

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

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
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
                "all_mapping_withFilter",
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

    bfragilis_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Bfragilis/coverage_filtered.csv"
    )
    ecoli_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Ecoli/coverage_filtered.csv"
    )

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield bfragilis_cov_fp, ecoli_cov_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    (
        bfragilis_cov_fp,
        ecoli_cov_fp,
        benchmarks_fp,
    ) = run_sunbeam

    # Check output
    assert os.path.exists(bfragilis_cov_fp)
    assert os.path.exists(ecoli_cov_fp)

    assert os.stat(bfragilis_cov_fp).st_size != 0
    assert os.stat(ecoli_cov_fp).st_size != 0