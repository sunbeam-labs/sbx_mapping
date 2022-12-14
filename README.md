<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_mapping

<!-- badges: start -->
[![Tests](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/tests.yml)
[![Super-Linter](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/linter.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_mapping/actions/workflows/linter.yml)
<!-- badges: end -->

A [Sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for mapping reads against a set of reference genomes.

## Installation

To install, activate your conda environment (using the name of your environment) and use `sunbeam extend`:

    conda activate <i>sunbeamX.X.X</i>
    sunbeam extend https://github.com/sunbeam-labs/sbx_mapping.git

## Usage

To generate a csv coverage file along with bam and bcf files, create a project, specify your references, and use the `all_mapping` target:

    sunbeam init --data_fp /path/to/reads/ /path/to/project/
    sunbeam config modify -i -f /path/to/project/sunbeam_config.yml -s 'sbx_mapping: {{genomes_fp: {/path/to/hosts/}}}'
    sunbeam run --profile /path/to/project/ all_mapping

N.B. For sunbeam versions <4 the last command will be something like `sunbeam run --configfile /path/to/project/sunbeam_config.yml all_mapping`.

## Configuration

  - genomes_fp: Is the filepath to your reference genomes (in fasta format)
  - samtools_opts: Are the options passed to samtools view

## Legacy Installation

For sunbeam versions <3 or if `sunbeam extend` isn't working, you can use `git` directly to install an extension:

    git clone https://github.com/sunbeam-labs/sbx_mapping.git extensions/sbx_mapping

and then include it in the config for any given project with:

    cat extensions/sbx_mapping/config.yml >> /path/to/project/sunbeam_config.yml