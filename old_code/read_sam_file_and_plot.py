import numpy as np
import matplotlib.pyplot as plt
import math


def get_bin_index_from_sam_line(line, quality_threshold):
    parts = line.split("\t")
    mapq_score = int(parts[4])
    if mapq_score >= quality_threshold:
        return int(round(int(parts[3]) / bin_size))
    else:
        return -1


def get_reads_count_per_position_from_file(genome_length, bin_size, my_mapping_file, quality_threshold):
    bin_size = bin_size * 1.0
    bins_count = int(math.ceil(genome_length / bin_size))
    coverage = np.array([0] * bins_count)
    infile = open(my_mapping_file, "r")
    line = infile.readline().rstrip()
    while line[0] == "@":
        line = infile.readline().rstrip()
    while line:
        bin_index =  get_bin_index_from_sam_line(line, quality_threshold)
        if bin_index != -1:
            coverage[bin_index] += 1
        line = infile.readline().rstrip()
    infile.close()
    return coverage


def save_fig_coverage(coverage, fig_name):
    plt.subplot(2, 1, 1)
    plt.title('whole chromosome ceoverage')
    plt.plot(coverage, 'o', markersize=3)
    plt.xticks([1000*x for x in range(7)], ['{} Mb'.format(x) for x in range(7)])
    plt.title('reads per Kb - chromosome')
    plt.ylabel('#reads')

    plt.subplot(2, 1, 2)
    plt.title('whole chromosome ceoverage')
    plt.plot(coverage[range(3000, 4000)], 'o', markersize=3)
    plt.xticks([x*200 for x in range(6)],  ['{} Kb'.format(3000+x*200) for x in range(6)])
    plt.title('reads per Kb - deletion')
    plt.ylabel('#reads')
    plt.tight_layout()
    plt.savefig(fig_name)
    plt.clf()


genome_length = 6472489
bin_size = 1000
quality_threshold = 0
my_mapping_file = #file name here#
coverage = get_reads_count_per_position_from_file(genome_length, bin_size, my_mapping_file, 0)
save_fig_coverage(coverage, 'all_reads.png')

coverage = get_reads_count_per_position_from_file(genome_length, bin_size, my_mapping_file, 2)
save_fig_coverage(coverage, 'th_reads.png')
