import sys

try:
    GenomeFiles
    GenomeSegments
except NameError:
    sys.stderr.write("sbx_mapping::INFO Collecting target genomes... ")
    if (
        Cfg["sbx_mapping"]["genomes_fp"] == Cfg["all"]["root"]
        or not Cfg["sbx_mapping"]["genomes_fp"]
    ):
        GenomeFiles = []
        GenomeSegments = {}
    else:
        GenomeFiles = [f for f in Cfg["sbx_mapping"]["genomes_fp"].glob("*.fasta")]
        GenomeSegments = {
            PurePath(g.name).stem: read_seq_ids(Cfg["sbx_mapping"]["genomes_fp"] / g)
            for g in GenomeFiles
        }
    sys.stderr.write("done.\n")
    sys.stderr.write(f"sbx_mapping::INFO Genome files found: {str(GenomeFiles)}\n")

TARGET_MAPPING_FILTER = [
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam.bai"),
        genome=GenomeSegments.keys(),
        sample=Samples.keys(),
    ),
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "coverage_filtered.csv"),
        genome=GenomeSegments.keys(),
    ),
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "numReads.csv"),
        genome=GenomeSegments.keys(),
    ),
    expand(str(MAPPING_FP / "{genome}" / "numReads.csv"), genome=GenomeSegments.keys()),
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "sliding_coverage.csv"),
        genome=GenomeSegments.keys(),
    ),
    expand(
        str(MAPPING_FP / "{genome}" / "sliding_coverage.csv"),
        genome=GenomeSegments.keys(),
    ),
]


rule get_sliding_coverage:
    input:
        str(MAPPING_FP / "{genome}" / "{sample}.bam"),
    output:
        str(MAPPING_FP / "intermediates" / "{genome}" / "{sample}_sliding_coverage.csv"),
    params:
        window_size=Cfg["sbx_mapping_withFilter"]["window_size"],
        sampling=Cfg["sbx_mapping_withFilter"]["sampling"],
    conda:
        "sbx_mapping_env.yml"
    script:
        "scripts/get_sliding_coverage.py"


def _sliding_coverage_csvs(w):
    pattern = str(
        MAPPING_FP / "intermediates" / w.genome / "{sample}_sliding_coverage.csv"
    )
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule summarize_sliding_coverage:
    input:
        _sliding_coverage_csvs,
    output:
        str(MAPPING_FP / "{genome}" / "sliding_coverage.csv"),
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"


rule all_mapping_withFilter:
    input:
        TARGET_MAPPING_FILTER,


rule filter_aln_quality:
    input:
        str(MAPPING_FP / "{genome}" / "{sample}.bam"),
    output:
        str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam"),
    params:
        alnLen=Cfg["sbx_mapping_withFilter"]["alnLen"],
        percIdentity=Cfg["sbx_mapping_withFilter"]["percIdentity"],
    conda:
        "sbx_mapping_env.yml"
    script:
        "scripts/filter_aln_quality.py"


# rule samtools_index_filtered:
#    input: str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam')
#    output: str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam.bai')
#    shell: "samtools index {input} {output}"


rule get_coverage_filtered:
    input:
        bam=str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam"),
        bai=str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam.bai"),
    output:
        str(MAPPING_FP / "filtered" / "intermediates" / "{genome}" / "{sample}.csv"),
    benchmark:
        BENCHMARK_FP / "get_coverage_filtered_{genome}_{sample}.tsv"
    log:
        LOG_FP / "get_coverage_filtered_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    script:
        "scripts/samtools_get_coverage.py"


def _sorted_filtered_csvs(w):
    pattern = str(MAPPING_FP / "filtered" / "intermediates" / w.genome / "{sample}.csv")
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule samtools_summarize_filtered_coverage:
    input:
        _sorted_filtered_csvs,
    output:
        str(MAPPING_FP / "filtered" / "{genome}" / "coverage_filtered.csv"),
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"


rule summarize_num_mapped_reads:
    input:
        str(MAPPING_FP / "{genome}" / "{sample}.bam"),
    output:
        str(MAPPING_FP / "intermediates" / "{genome}" / "{sample}_numReads.csv"),
    shell:
        """
        samtools idxstats {input} | (sed 's/^/{wildcards.sample}\t/') > {output}
        """


def _numReads(w):
    pattern = str(MAPPING_FP / "intermediates" / w.genome / "{sample}_numReads.csv")
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule summarize_num_reads:
    input:
        _numReads,
    output:
        str(MAPPING_FP / "{genome}" / "numReads.csv"),
    shell:
        "(cat {input}) > {output}"
