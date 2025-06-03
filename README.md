<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_mapping

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![Release](https://img.shields.io/github/release/sunbeam-labs/sbx_mapping.svg?style=flat)](https://github.com/sunbeam-labs/sbx_mapping/releases/latest)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_mapping)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_mapping/)
<!-- badges: end -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for mapping reads against a set of reference genomes.

## Config

  - genomes_fp: Is the filepath to your reference genomes (in fasta format) **all references must have .fasta OR .fa OR .fna extension!!**
  - samtools_opts: Are the options passed to samtools view
  - alnLen: Is the minimum alignment length
  - percIdentity: Is the minimum percent identity
  - window_size: Is the size of the sliding window for coverage calculations
  - sampling: Is the sampling rate for coverage calculations

## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).