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


def 

import os
reads_fasta = '/prodslow/testing/ariel/genomagic_qa/SGU-29/nam_6k_design.fasta'
ref_fasta = '/prodslow/testing/ariel/genomagic_qa/SGU-29/williams82__ver100.fasta'
verify_indexing_of_fasta(ref_fasta)

#csv_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/out.csv'
#filtered_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/filtered.csv'
#sam_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/a.sam'
#f = open(sam_file, 'r')
#csv = open(csv_file,'w')
#csv_header = "marker_id,genome_id,chr_id,chr_start,chr_end,chr_strand,alg_sig,alg_map_quality,seqA,seqB,isSeqA_informative,isSeqB_informative"
#csv.write("{}\n".format(csv_header))

#filtered = open(filtered_file,'w')
#allele1_aln = get_next_sam_line_as_list(f)
#while allele1_aln is not None:
#    allele2_aln = get_next_sam_line_as_list(f)
#    probe_id1 = allele1_aln[0]
#    probe_id2 = allele2_aln[0]
#    assert probe_id1[-2:] == '_1'
#    assert probe_id2[-2:] == '_2'
#    assert probe_id1[:-2] == probe_id2[:-2]
#    handle_marker_alignment(allele1_aln, allele2_aln, csv, filtered)
#    allele1_aln = get_next_sam_line_as_list(f)
#f.close()
#csv.close()
#filtered.close()


