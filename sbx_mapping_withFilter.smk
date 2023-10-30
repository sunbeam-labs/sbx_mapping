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


rule samtools_get_sliding_coverage:
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
        "scripts/samtools_get_sliding_coverage.py"


def _sliding_coverage_csvs(w):
    pattern = str(
        MAPPING_FP / "intermediates" / w.genome / "{sample}_sliding_coverage.csv"
    )
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule samtools_summarize_sliding_coverage:
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
    run:
        print(input)
        filter_bam_alignments(input[0], output[0], params.percIdentity, params.alnLen)


# rule samtools_index_filtered:
#    input: str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam')
#    output: str(MAPPING_FP/'filtered'/'{genome}'/'{sample}.bam.bai')
#    shell: "samtools index {input} {output}"


rule samtools_get_coverage_filtered:
    input:
        in_bam=str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam"),
        in_bai=str(MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam.bai"),
    output:
        str(MAPPING_FP / "filtered" / "intermediates" / "{genome}" / "{sample}.csv"),
    run:
        samtools.get_coverage_stats(
            wildcards.genome, input.in_bam, wildcards.sample, output[0]
        )


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


rule samtools_summarize_num_mapped_reads:
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


rule samtools_summarize_numReads:
    input:
        _numReads,
    output:
        str(MAPPING_FP / "{genome}" / "numReads.csv"),
    shell:
        "(cat {input}) > {output}"
