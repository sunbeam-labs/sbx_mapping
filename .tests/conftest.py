import os
import pytest
import shutil
import subprocess as sp
import sys


@pytest.fixture
def setup(tmpdir):
    reads_fp = os.path.abspath(".tests/data/reads/")
    genomes_fp = os.path.abspath(".tests/data/ref/")

    project_dir = os.path.join(tmpdir, "project/")

    profile = os.environ.get("TEST_PROFILE", "")
    if profile:
        profile_list = ["--profile", profile]
    else:
        profile_list = []

    sp.check_output(
        ["sunbeam", "init", "--data_fp", reads_fp, project_dir] + profile_list
    )

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

    yield tmpdir, project_dir

    shutil.rmtree(tmpdir)


@pytest.fixture
def run_sunbeam(setup):
    tmpdir, project_dir = setup

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
                tmpdir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
        shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")
        sys.exit(e)

    shutil.copytree(os.path.join(output_fp, "logs/"), "logs/")
    shutil.copytree(os.path.join(project_dir, "stats/"), "stats/")

    bfragilis_sliding_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Bfragilis/sliding_coverage.csv"
    )
    ecoli_sliding_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Ecoli/sliding_coverage.csv"
    )
    bfragilis_filtered_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Bfragilis/coverage_filtered.csv"
    )
    ecoli_filtered_cov_fp = os.path.join(
        output_fp, "mapping/filtered/Ecoli/coverage_filtered.csv"
    )
    bfragilis_num_reads_fp = os.path.join(
        output_fp, "mapping/filtered/Bfragilis/numReads.tsv"
    )
    ecoli_num_reads_fp = os.path.join(output_fp, "mapping/filtered/Ecoli/numReads.tsv")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield bfragilis_sliding_cov_fp, ecoli_sliding_cov_fp, bfragilis_filtered_cov_fp, ecoli_filtered_cov_fp, bfragilis_num_reads_fp, ecoli_num_reads_fp, benchmarks_fp
