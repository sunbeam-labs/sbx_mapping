# -*- mode: Snakemake -*-

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

TARGET_MAPPING = [
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "coverage_filtered.csv"),
        genome=GenomeSegments.keys(),
    ),
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "numReads.csv"),
        genome=GenomeSegments.keys(),
    ),
    expand(
        str(MAPPING_FP / "filtered" / "{genome}" / "sliding_coverage.csv"),
        genome=GenomeSegments.keys(),
    ),
]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


ruleorder: build_host_index > build_genome_index


localrules:
    all_mapping,


rule all_mapping:
    input:
        TARGET_MAPPING,


### Prepping BAM files ###


rule build_genome_index:
    input:
        Cfg["sbx_mapping"]["genomes_fp"] / "{genome}.fasta",
    output:
        [
            Cfg["sbx_mapping"]["genomes_fp"] / ("{genome}.fasta." + ext)
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
    benchmark:
        BENCHMARK_FP / "build_genome_index_{genome}.tsv"
    log:
        LOG_FP / "build_genome_index_{genome}.log",
    conda:
        "sbx_mapping_env.yml"
    shell:
        "cd {Cfg[sbx_mapping][genomes_fp]} && bwa index {input} 2>&1 | tee {log}"


rule align_to_genome:
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        index=Cfg["sbx_mapping"]["genomes_fp"] / "{genome}.fasta.amb",
    output:
        temp(MAPPING_FP / "intermediates" / "{genome}" / "{sample}.sam"),
    benchmark:
        BENCHMARK_FP / "align_to_genome_{genome}_{sample}.tsv"
    log:
        LOG_FP / "align_to_genome_{genome}_{sample}.log",
    params:
        index_fp=Cfg["sbx_mapping"]["genomes_fp"],
    threads: 4
    conda:
        "sbx_mapping_env.yml"
    shell:
        """
        bwa mem -M -t {threads} \
        {params.index_fp}/{wildcards.genome}.fasta \
        {input.reads} -o {output} \
        2>&1 | tee {log}
        """


rule samtools_convert:
    input:
        MAPPING_FP / "intermediates" / "{genome}" / "{sample}.sam",
    output:
        MAPPING_FP / "{genome}" / "{sample}.bam",
    benchmark:
        BENCHMARK_FP / "samtools_convert_{genome}_{sample}.tsv"
    log:
        view_log=LOG_FP / "samtools_convert_view_{genome}_{sample}.log",
        sort_log=LOG_FP / "samtools_convert_sort_{genome}_{sample}.log",
    threads: 4
    conda:
        "sbx_mapping_env.yml"
    shell:
        """
        samtools view -@ {threads} -b {Cfg[sbx_mapping][samtools_opts]} {input} 2>&1 | tee {log.view_log} | \
        samtools sort -@ {threads} -o {output} -O bam 2>&1 | tee {log.sort_log}
        """


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


rule samtools_index:
    input:
        MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam.bai",
    benchmark:
        BENCHMARK_FP / "samtools_getindex_{genome}_{sample}.tsv"
    log:
        LOG_FP / "samtools_index_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    shell:
        "samtools index {input} {output} 2>&1 | tee {log}"


### Calculate sliding window coverage stats ###


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


### Get filtered coverage stats ###

        
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


### Get number of mapped reads ###


rule summarize_num_mapped_reads:
    input:
        str(MAPPING_FP / "{genome}" / "{sample}.bam"),
    output:
        str(MAPPING_FP / "intermediates" / "{genome}" / "{sample}_numReads.csv"),
    conda:
        "sbx_mapping_env.yml"
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




rule samtools_mpileup:
    input:
        bam=MAPPING_FP / "{genome}" / "{sample}.bam",
        genome=Cfg["sbx_mapping"]["genomes_fp"] / "{genome}.fasta",
    output:
        mpileup_out=temp(MAPPING_FP / "{genome}" / "{sample}.bcf"),
        call_out=MAPPING_FP / "{genome}" / "{sample}.raw.bcf",
    benchmark:
        BENCHMARK_FP / "samtools_mpileup_{genome}_{sample}.tsv"
    log:
        mpileup_log=LOG_FP / "samtools_mpileup_mpileup_{genome}_{sample}.log",
        call_log=LOG_FP / "samtools_mpileup_call_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    shell:
        """
        bcftools mpileup -f {input.genome} {input.bam} -o {output.mpileup_out} 2>&1 | tee {log.mpileup_log} && \
        bcftools call -Ob -v -c -o {output.call_out} {output.mpileup_out} 2>&1 | tee {log.call_log}
        """
