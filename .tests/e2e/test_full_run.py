import os
import shutil
import subprocess as sp
import tempfile
import unittest


class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.reads_fp = os.path.abspath(".tests/data/reads/")
        self.genomes_fp = os.path.abspath(".tests/data/hosts/")

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output(
            ["sunbeam", "init", "--data_fp", self.reads_fp, self.project_dir]
        )

        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        config_str = f"sbx_mapping: {{genomes_fp: {self.genomes_fp}}}"

        sp.check_output(
            [
                "sunbeam",
                "config",
                "modify",
                "-i",
                "-s",
                f"{config_str}",
                f"{self.config_fp}",
            ]
        )

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")
        # shutil.copytree(".tests/data/sunbeam_output", self.output_fp)

        self.human_genome_fp = os.path.join(self.output_fp, "mapping/human/")
        self.human_copy_genome_fp = os.path.join(self.output_fp, "mapping/human_copy/")
        self.phix174_genome_fp = os.path.join(self.output_fp, "mapping/phix174/")

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_full_run(self):
        # Run the test job.
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                self.project_dir,
                "all_mapping",
                "--directory",
                self.temp_dir,
            ]
        )

        # Check output
        self.assertTrue(os.listdir(self.human_genome_fp) == [])
        self.assertTrue(os.listdir(self.human_copy_genome_fp) == [])
        self.assertTrue(os.listdir(self.phix174_genome_fp) == [])
