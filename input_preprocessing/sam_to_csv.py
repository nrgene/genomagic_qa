supplementary_alignment_flag = 2048
max_dist_alleles_alignment = 30
min_score_threshold = 20

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
    same_chromosome = allele1_aln[2] == allele2_aln[2]
    seq_proximity = abs(int(allele1_aln[3]) - int(allele1_aln[3])) <= max_dist_alleles_alignment
    my_score = min(int(allele1_aln[4]), int(allele2_aln[4]))
    scores_are_high = my_score >= min_score_threshold
    same_seq = ''





csv_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/out.csv'
filtered_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/filtered.csv'
sam_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/a.sam'
f = open(sam_file, 'r')
csv = open(csv_file,'w')
filtered = open(filtered_file,'w')
allele1_aln = get_next_sam_line_as_list(f)
while allele1_aln is not None:
    allele2_aln = get_next_sam_line_as_list(f)
    probe_id1 = allele1_aln[0]
    probe_id2 = allele2_aln[0]
    assert probe_id1[-2:] == '_1'
    assert probe_id2[-2:] == '_2'
    assert probe_id1[:-2] == probe_id2[:-2]
    handle_marker_alignment(allele1_aln, allele2_aln, csv_file, filtered_file)
    allele1_aln = get_next_sam_line_as_list(f)
f.close()
csv.close()
filtered.close()


