<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_mapping

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml)
[![Release](https://img.shields.io/github/release/sunbeam-labs/sbx_mapping.svg?style=flat)](https://github.com/sunbeam-labs/sbx_mapping/releases/latest)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_mapping)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_mapping/)
<!-- badges: end -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for mapping reads against a set of reference genomes.

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_mapping.git

## Usage

To generate coverage reports, create a project, specify your references, and use the `all_mapping` target:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    sunbeam config modify -i -f /path/to/project/sunbeam_config.yml -s 'sbx_mapping: {{genomes_fp: {/path/to/uncompressed/fasta/reference/genomes/}}}'
    sunbeam run --profile /path/to/project/ all_mapping

N.B. For sunbeam versions <4 the last command will be something like `sunbeam run --configfile /path/to/project/sunbeam_config.yml all_mapping`.

## Configuration

  - genomes_fp: Is the filepath to your reference genomes (in fasta format) **all references must have .fasta OR .fa OR .fna extension!!**
  - samtools_opts: Are the options passed to samtools view
  - alnLen: Is the minimum alignment length
  - percIdentity: Is the minimum percent identity
  - window_size: Is the size of the sliding window for coverage calculations
  - sampling: Is the sampling rate for coverage calculations

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_mapping.git extensions/sbx_mapping

and then include it in the config for any given project with:

    cat extensions/sbx_mapping/config.yml >> /path/to/project/sunbeam_config.yml
