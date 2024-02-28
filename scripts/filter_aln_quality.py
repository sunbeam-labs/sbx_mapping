# from scripts.filter_aln_quality_f import filter_bam_alignments
import pysam
import re


def filter_bam_alignments(bam_in_fp, bam_out_fp, percIdentity, alnLen):
    """
    Filter the alignments in the bam file with the defined percent identity
    and alingment length thresholds.
    bam_in_fp:BAM file name to filer
    bam_out_fp:BAM file name to output
    percIdentity:percent identity threshold (out of 1)
    alnLen:alignment length threshold
    """
    f_in = pysam.AlignmentFile(bam_in_fp)
    with pysam.AlignmentFile(bam_out_fp, "wb", template=f_in) as out_file:
        for item in f_in:
            if item.has_tag("MD"):
                mdstr = item.get_tag("MD")
                mdSub = re.sub(r"([\\^]*[ACGT]+)[0]*", " \\1 ", mdstr)
                mdSplit = re.split("[ ]+", mdSub)
                nums = [int(i) for i in mdSplit if i.isdigit()]
                letters = [i for i in mdSplit if not i.isdigit()]
                letters = re.sub("[^ATCG]", "", "".join(letters))

                alnLen_seq = sum(nums) + len(letters)
                percIdentity_seq = sum(nums) / alnLen_seq

                if alnLen_seq > alnLen and percIdentity_seq > percIdentity:
                    out_file.write(item)
    f_in.close()


filter_bam_alignments(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.params.percIdentity,
    snakemake.params.alnLen,
)
