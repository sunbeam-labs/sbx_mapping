# -*- mode: Snakemake -*-

import sys

sys.stderr.write("Collecting target genomes... ")
if Cfg["sbx_mapping"]["genomes_fp"] == Cfg["all"]["root"]:
    GenomeFiles = []
    GenomeSegments = {}
else:
    GenomeFiles = [f for f in Cfg["sbx_mapping"]["genomes_fp"].glob("*.fasta")]
    GenomeSegments = {
        PurePath(g.name).stem: read_seq_ids(Cfg["sbx_mapping"]["genomes_fp"] / g)
        for g in GenomeFiles
    }
sys.stderr.write("done.\n")

TARGET_MAPPING = [
    expand(
        MAPPING_FP / "{genome}" / "{sample}.bam.bai",
        genome=GenomeSegments.keys(),
        sample=Samples.keys(),
    ),
    expand(
        MAPPING_FP / "{genome}" / "{sample}.raw.bcf",
        genome=GenomeSegments.keys(),
        sample=Samples.keys(),
    ),
    expand(MAPPING_FP / "{genome}" / "coverage.csv", genome=GenomeSegments.keys()),
]


ruleorder: build_host_index > build_genome_index


localrules:
    all_mapping,


rule all_mapping:
    input:
        TARGET_MAPPING,


rule build_genome_index:
    input:
        Cfg["sbx_mapping"]["genomes_fp"] / "{genome}.fasta",
    output:
        Cfg["sbx_mapping"]["genomes_fp"] / "{genome}.fasta.amb",
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


def _sorted_csvs(w):
    pattern = MAPPING_FP / "intermediates" / w.genome / "{sample}.csv"
    paths = sorted(expand(pattern, sample=Samples.keys()))
    return paths


rule samtools_summarize_coverage:
    input:
        _sorted_csvs,
    output:
        MAPPING_FP / "{genome}" / "coverage.csv",
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output}"


rule samtools_get_coverage:
    input:
        MAPPING_FP / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP / "intermediates" / "{genome}" / "{sample}.csv",
    benchmark:
        BENCHMARK_FP / "samtools_get_coverage_{genome}_{sample}.tsv"
    log:
        LOG_FP / "samtools_get_coverage_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    script:
        "scripts/samtools_get_coverage.py"


rule samtools_index:
    input:
        MAPPING_FP / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP / "{genome}" / "{sample}.bam.bai",
    benchmark:
        BENCHMARK_FP / "samtools_getindex_{genome}_{sample}.tsv"
    log:
        LOG_FP / "samtools_index_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    shell:
        "samtools index {input} {output} 2>&1 | tee {log}"


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
        mpileup_log=LOG_FP / "samtools_mpileup_{genome}_{sample}.log",
        call_log=LOG_FP / "samtools_mpileup_{genome}_{sample}.log",
    conda:
        "sbx_mapping_env.yml"
    shell:
        """
        bcftools mpileup -f {input.genome} {input.bam} -o {output.mpileup_out} 2>&1 | tee {log.mpileup_log} && \
        bcftools call -Ob -v -c -o {output.call_out} {output.mpileup_out} 2>&1 | tee {log.call_log}
        """
