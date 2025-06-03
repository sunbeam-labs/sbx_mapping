from pathlib import PurePath
from sunbeam.bfx.parse import parse_fasta


try:
    SBX_MAPPING_VERSION = get_ext_version("sbx_mapping")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_MAPPING_VERSION = "0.0.0"

if (
    Cfg["qc"]["host_fp"] == Cfg["sbx_mapping"]["genomes_fp"]
    and Cfg["qc"]["host_fp"] != Cfg["all"]["root"]
):
    raise ValueError("sbx_mapping::ERROR: Host and target genomes cannot be the same")


def read_seq_ids(fasta_fp: str) -> List[Tuple[str, str]]:
    """
    Return the sequence identifiers for a given fasta filepath.
    """
    with open(fasta_fp) as f:
        return list(parse_fasta(f))


HOST_FILE_EXT = ".fasta"
sys.stderr.write("sbx_mapping::INFO Collecting target genomes... ")
if (
    Cfg["sbx_mapping"]["genomes_fp"] == Cfg["all"]["root"]
    or not Cfg["sbx_mapping"]["genomes_fp"]
):
    GenomeFiles = []
    GenomeSegments = {}
else:
    GenomeFiles = [f for f in Cfg["sbx_mapping"]["genomes_fp"].glob("*.fasta")]
    if not GenomeFiles:
        GenomeFiles = [f for f in Cfg["sbx_mapping"]["genomes_fp"].glob("*.fa")]
        HOST_FILE_EXT = ".fa"
    if not GenomeFiles:
        GenomeFiles = [f for f in Cfg["sbx_mapping"]["genomes_fp"].glob("*.fna")]
        HOST_FILE_EXT = ".fna"
    GenomeSegments = {
        PurePath(g.name).stem: read_seq_ids(Cfg["sbx_mapping"]["genomes_fp"] / g)
        for g in GenomeFiles
    }
    GenomeFiles = {PurePath(g.name).stem: g for g in GenomeFiles}
sys.stderr.write("done.\n")
sys.stderr.write(f"sbx_mapping::INFO Genome files found: {str(GenomeFiles)}\n")


localrules:
    all_mapping,


rule all_mapping:
    input:
        expand(
            MAPPING_FP / "filtered" / "{genome}" / "coverage_filtered.csv",
            genome=GenomeSegments.keys(),
        ),
        expand(
            MAPPING_FP / "filtered" / "{genome}" / "numReads.tsv",
            genome=GenomeSegments.keys(),
        ),
        expand(
            MAPPING_FP / "filtered" / "{genome}" / "sliding_coverage.csv",
            genome=GenomeSegments.keys(),
        ),


rule build_genome_index:
    input:
        Cfg["sbx_mapping"]["genomes_fp"] / ("{genome}" + HOST_FILE_EXT),
    output:
        [
            Cfg["sbx_mapping"]["genomes_fp"] / ("{genome}" + HOST_FILE_EXT + "." + ext)
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
    benchmark:
        BENCHMARK_FP / "build_genome_index_{genome}.tsv"
    log:
        LOG_FP / "build_genome_index_{genome}.log",
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    shell:
        """
        cd $(dirname {input})
        bwa index {input} > {log} 2>&1
        """


rule align_to_genome:
    input:
        *rules.build_genome_index.output,
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        temp(MAPPING_FP / "intermediates" / "{genome}" / "{sample}.sam"),
    benchmark:
        BENCHMARK_FP / "align_to_genome_{genome}_{sample}.tsv"
    log:
        LOG_FP / "align_to_genome_{genome}_{sample}.log",
    params:
        index_fp=Cfg["sbx_mapping"]["genomes_fp"],
        host_file_ext=HOST_FILE_EXT,
    threads: 4
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    shell:
        """
        bwa mem -M -t {threads} \
        {params.index_fp}/{wildcards.genome}{params.host_file_ext} \
        {input.reads} -o {output} \
        > {log} 2>&1
        """


rule samtools_convert:
    input:
        MAPPING_FP / "intermediates" / "{genome}" / "{sample}.sam",
    output:
        MAPPING_FP / "{genome}" / "{sample}.bam",
    benchmark:
        BENCHMARK_FP / "samtools_convert_{genome}_{sample}.tsv"
    log:
        LOG_FP / "samtools_convert_{genome}_{sample}.log",
    threads: 4
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    shell:
        """
        (
            samtools view -@ {threads} -b {Cfg[sbx_mapping][samtools_opts]} {input} | \
            samtools sort -@ {threads} -o {output} -O bam
        ) > {log} 2>&1
        """


rule filter_aln_quality:
    input:
        MAPPING_FP / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam",
    params:
        alnLen=Cfg["sbx_mapping"]["alnLen"],
        percIdentity=Cfg["sbx_mapping"]["percIdentity"],
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
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
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    shell:
        "samtools index {input} {output} > {log} 2>&1"


rule get_sliding_coverage:
    input:
        MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP
        / "filtered"
        / "intermediates"
        / "{genome}"
        / "{sample}_sliding_coverage.csv",
    params:
        window_size=Cfg["sbx_mapping"]["window_size"],
        sampling=Cfg["sbx_mapping"]["sampling"],
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    script:
        "scripts/get_sliding_coverage.py"


rule summarize_sliding_coverage:
    input:
        sorted(
            expand(
                MAPPING_FP
                / "filtered"
                / "intermediates"
                / "{{genome}}"
                / "{sample}_sliding_coverage.csv",
                sample=Samples.keys(),
            )
        ),
    output:
        MAPPING_FP / "filtered" / "{genome}" / "sliding_coverage.csv",
    benchmark:
        BENCHMARK_FP / "summarize_sliding_coverage_{genome}.tsv"
    log:
        LOG_FP / "summarize_sliding_coverage_{genome}.log",
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output} 2> {log}"


rule get_coverage_filtered:
    input:
        bam=MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam",
        bai=MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam.bai",
    output:
        MAPPING_FP / "filtered" / "coverage" / "{genome}" / "{sample}.csv",
    benchmark:
        BENCHMARK_FP / "get_coverage_filtered_{genome}_{sample}.tsv"
    log:
        LOG_FP / "get_coverage_filtered_{genome}_{sample}.log",
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    script:
        "scripts/samtools_get_coverage.py"


rule samtools_summarize_filtered_coverage:
    input:
        sorted(
            expand(
                MAPPING_FP / "filtered" / "coverage" / "{{genome}}" / "{sample}.csv",
                sample=Samples.keys(),
            )
        ),
    output:
        MAPPING_FP / "filtered" / "{genome}" / "coverage_filtered.csv",
    benchmark:
        BENCHMARK_FP / "samtools_summarize_filtered_coverage_{genome}.tsv"
    log:
        LOG_FP / "samtools_summarize_filtered_coverage_{genome}.log",
    shell:
        "(head -n 1 {input[0]}; tail -q -n +2 {input}) > {output} 2> {log}"


rule summarize_num_mapped_reads:
    input:
        MAPPING_FP / "filtered" / "{genome}" / "{sample}.bam",
    output:
        MAPPING_FP / "filtered" / "intermediates" / "{genome}" / "{sample}_numReads.tsv",
    conda:
        "envs/sbx_mapping_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_mapping:{SBX_MAPPING_VERSION}"
    shell:
        """
        (
            samtools idxstats {input} | \
            (sed 's/^/{wildcards.sample}\t/') \
            > {output}
        ) > {log} 2>&1
        """


rule summarize_num_reads:
    input:
        sorted(
            expand(
                MAPPING_FP
                / "filtered"
                / "intermediates"
                / "{{genome}}"
                / "{sample}_numReads.tsv",
                sample=Samples.keys(),
            )
        ),
    output:
        MAPPING_FP / "filtered" / "{genome}" / "numReads.tsv",
    shell:
        "(cat {input}) > {output}"
