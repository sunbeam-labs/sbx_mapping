import sys
from sunbeamlib import samtools

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    samtools.get_coverage_stats(
        snakemake.wildcards.genome,
        snakemake.input[0],
        snakemake.wildcards.sample,
        snakemake.output[0],
    )
