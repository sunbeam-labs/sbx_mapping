import csv
import numpy
import subprocess


def sliding_window_coverage(genome, bamfile, sample, output_fp, N, sampling):
    args = ["samtools", "depth", "-aa", bamfile]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, universal_newlines=True)
    # Organize into a list of depths for each segment, streaming in text
    reader = csv.reader(p.stdout, delimiter="\t")
    data = {}
    for row in reader:
        if not data.get(row[0]):
            data[row[0]] = []
        data[row[0]].append(int(row[2]))

    fields = ["Genome", "Segment", "Sample", "Location", "Average"]
    with open(output_fp, "w") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for segment in data.keys():
            if len(data[segment]) > sampling:
                moving_avg = numpy.convolve(
                    data[segment], numpy.ones((N,)) / N, mode="full"
                )
                for i, x in enumerate(moving_avg):
                    if i % sampling == 0:
                        writer.writerow([genome, segment, sample, i, x])


sliding_window_coverage(
    snakemake.wildcards.genome,
    snakemake.input[0],
    snakemake.wildcards.sample,
    snakemake.output[0],
    snakemake.params.window_size,
    snakemake.params.sampling,
)
