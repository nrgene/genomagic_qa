import os

supplementary_alignment_flag = 2048
max_dist_alleles_alignment = 30
min_score_threshold = 20
min_sam_columns_num = 13

def get_next_sam_line_as_list(f):
    while True:
        next_line = f.readline().rstrip()
        if not next_line:
            return None
        if next_line[0] != '@':
            arr = next_line.split('\t')
            curr_flag = int(arr[1])
            if curr_flag < supplementary_alignment_flag:
                return arr

def handle_marker_alignment(allele1_aln, allele2_aln, out_include, out_exclude):
    assert allele1_aln[1] == "0" or allele1_aln[1] == "16" or allele1_aln[1] == "4", allele1_aln[1]
    assert allele2_aln[1] == "0" or allele2_aln[1] == "16" or allele2_aln[1] == "4", allele2_aln[1]
    same_chromosome = allele1_aln[2] == allele2_aln[2]
    seq_proximity = abs(int(allele1_aln[3]) - int(allele1_aln[3])) <= max_dist_alleles_alignment
    my_score = min(int(allele1_aln[4]), int(allele2_aln[4]))
    scores_are_high = my_score >= min_score_threshold
    if same_chromosome and seq_proximity and scores_are_high:
        out_include.write("")
        print('yes')
    else:
        print('no')
#    same_flag = ''



def chech_if_ref_has_index(ref_fasta):
    ref_fasta_suffix = ref_fasta.split('.')[-1]
    assert ref_fasta_suffix == "fasta" or ref_fasta_suffix == "fa", "bad format ({}), expecting fa/fasta".format(
        ref_fasta_suffix)
    index_file_name = "{}.bwt".format(ref_fasta)
    has_index = os.path.exists(index_file_name)
    return has_index


def verify_indexing_of_fasta(ref_fasta):
    has_index = chech_if_ref_has_index(ref_fasta)
    if has_index:
        print("found ref index")
    else:
        print("missing ref index. Indexing..")
        os.system("bwa index {}".format(ref_fasta))
    assert chech_if_ref_has_index(ref_fasta)

def bwa_align(reads_fasta, ref_fasta, reads_sam, threads_num):
    verify_indexing_of_fasta(ref_fasta)
    print('running BWA on genome {}'.format(ref_fasta))
    bwa_command = 'bwa mem -t {} {} {} | grep -v @ | awk \'$2<2048\'> {}'.format(threads_num, ref_fasta, reads_fasta, reads_sam)
    os.system(bwa_command)

def verify_seq_names(probe1_id, probe2_id):
    assert probe1_id[-2:] == "_1", probe1_id[-2:]
    assert probe2_id[-2:] == "_2", probe2_id[-2:]
    assert probe1_id[:-2] == probe2_id[:-2]


def verify_flags(curr_flag1, curr_flag2):
    assert curr_flag1 == "0" or curr_flag1 == "16" or curr_flag1 == "4", curr_flag1
    assert curr_flag2 == "0" or curr_flag2 == "16" or curr_flag2 == "4", curr_flag2


def handle_two_sam_lines(allele1_line, allele2_line, genome_name, out_include, out_exclude):
    al1_arr = allele1_line.split('\t')
    al2_arr = allele2_line.split('\t')
    assert len(al1_arr) >= min_sam_columns_num
    assert len(al2_arr) >= min_sam_columns_num
    probe1_id = al1_arr[0]
    probe2_id = al2_arr[0]
    verify_seq_names(probe1_id, probe2_id)
    flag1 = al1_arr[1]
    flag2 = al2_arr[1]
    verify_flags(flag1, flag2)
    chr1 = al1_arr[2]
    chr2 = al2_arr[2]
    same_chromosome = chr1 == chr2
    pos1 = int(al1_arr[3])
    pos2 = int(al2_arr[3])
    seq_proximity = abs(pos1 - pos2) <= max_dist_alleles_alignment
    both_scores_are_high = min(int(al1_arr[4]), int(al2_arr[4])) >= min_score_threshold
    seq1 = al1_arr[8]
    seq2 = al2_arr[8]
    write_to_csv = same_chromosome and seq_proximity and both_scores_are_high
    if write_to_csv:
        # marker_id,genome_id
        out_include.write("{},{}".format(probe1_id[:-2], genome_name))
        # chr_id,chr_start,chr_end,chr_strand
        out_include.write(",{},{},{},+".format(chr1, pos1, pos1 + len(seq1)))
        # alg_sig,alg_map_quality,seqA,seqB,isSeqA_informative,isSeqB_informative
        out_include.write(",{},{},{},{},1,1\n".format(al1_arr[5], al1_arr[4], seq1, seq2))
    else:
        out_exclude.write(">{}\n{}\n>{}\n{}\n".format(probe1_id, seq1, probe2_id, seq2))
    return write_to_csv


def split_probes_to_csv_and_fasta(reads_sam, csv, exclude_fasta, genome_name):
    included_count = 0
    excluded_count = 0
    excluded = open(exclude_fasta,'w')
    f = open(reads_sam, 'r')
    line1 = f.readline().rstrip()
    while line1:
        line2 = f.readline().rstrip()
        is_included = handle_two_sam_lines(line1, line2, genome_name,csv, excluded)
        if is_included:
            included_count += 1
        else:
            excluded_count += 1
        line1 = f.readline().rstrip()
    f.close()
    excluded.close()
    print("{}: {} probes were added to csv, {} were written to fasta".format(genome_name, included_count, excluded_count))

# we get the fasta input, align, and each seq is splitted either to the csv or to the filtered_fasta file
def align_and_write_to_csv(input_fasta, ref_fasta, ref_name, csv_file, filtered_fasta, threads_num):
    reads_sam = '{}.sam'.format(input_fasta)
    bwa_align(input_fasta, ref_fasta, reads_sam, threads_num)
    split_probes_to_csv_and_fasta(reads_sam, csv_file, filtered_fasta, ref_name)


def align_to_multiple_genomes_and_write_csv(csv_file, ref_file_list, ref_name_list, input_fasta, threads_num):
    genomes_count = len(ref_file_list)
    assert len(ref_name_list) == genomes_count
    csv = open(csv_file, 'w')
    csv_header = "marker_id,genome_id,chr_id,chr_start,chr_end,chr_strand,alg_sig,alg_map_quality,seqA,seqB,isSeqA_informative,isSeqB_informative"
    csv.write("{}\n".format(csv_header))
    next_fasta_file = input_fasta
    for i in range(genomes_count):
        curr_ref_fasta = ref_file_list[i]
        curr_ref_name = ref_name_list[i]
        excluded_fasta = "{}.excluded.{}.fasta".format(input_fasta, i)
        align_and_write_to_csv(next_fasta_file, curr_ref_fasta, curr_ref_name, csv, excluded_fasta, threads_num)
        next_fasta_file = excluded_fasta
    csv.close()

