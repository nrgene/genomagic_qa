import os
import sys
working_dir_path = os.getcwd()
genomagic_qa_repo_path = '/'.join(working_dir_path.split('/')[:-1])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_iterator as vi
vcf_info_columns_num = 9


def get_position_from_vcf_splitted_line(parts):
    chr1 = int(parts[0])
    start1 = int(parts[1])
    assert parts[7][:7] == 'HS;END='
    end1 = int(parts[7][7:])
    return [chr1, start1, end1]


def get_curr_value(index_in_parents_vcf, parents_data_parts, nahap_value, seg_len):
    print(parents_data_parts[4])
    if index_in_parents_vcf == -1:
        return '{}:0,0,0:{}'.format(nahap_value, seg_len)
    else:
        new_val = parents_data_parts[vcf_info_columns_num + index_in_parents_vcf].split(':')
        assert len(new_val) == 3
        if new_val[0] == nahap_value:
            assert new_val[1] == 'in {} 0,0,0', '{} should be 0,0,0'.format(new_val[0], new_val[1])
            return '{}:0,0,0:{}'.format(nahap_value, seg_len)
        else:
            return '{}:{}:{}'.format(new_val[0], new_val[1], seg_len)


def get_line_info(chr_name, start, end, haps_list):
    curr_str = '{}\t{}\t.\tN\t{}\t1\tPASS\tHS;END={}\tGT:CO1:LN'.format(chr_name, start, haps_list, end)
    return curr_str


def get_line_data(progeny_data_parts, index_from_progeny_to_parents, parents_data_parts, seg_len):
    genotype_data_line=''
    nahap_ind = parents_data_parts[4].count(',') + 1
    nahap_value = '{}|{}'.format(nahap_ind, nahap_ind)
    for g in progeny_data_parts[vcf_info_columns_num:]:
        nv = get_new_value(g, index_from_progeny_to_parents, nahap_value, parents_data_parts, seg_len)
        genotype_data_line = '{}\t{}'.format(genotype_data_line, nv)
    return genotype_data_line


def get_new_value(old_value, index_from_progeny_to_parents, nahap_value, parents_data_parts, seg_len):
    curr_genotype = old_value.split('|')[0]
    assert curr_genotype in index_from_progeny_to_parents
    index_in_parents_vcf = index_from_progeny_to_parents[curr_genotype]
    if index_in_parents_vcf == -1:
        new_val = "{}:0,0,0:{}".format(nahap_value, seg_len)
    else:
        assert index_in_parents_vcf >= 0, index_in_parents_vcf
        parent_val = parents_data_parts[vcf_info_columns_num + index_in_parents_vcf].split(':')
        if parent_val[0] == nahap_value:
            assert parent_val[1] == '0,0,0', parent_val
            new_val = "{}:0,0,0:{}".format(nahap_value, seg_len)
        else:
            new_val = "{}:{}:{}".format(parent_val[0], parent_val[1], seg_len)
    return new_val


def get_mapping_from_progeny_to_parents(haps_list_of_progeny, progeny_samples, parents_samples):
    coloring_samples_count = haps_list_of_progeny.count(',')
    mapping_from_progeny_to_paretns = {}
    for i in range(coloring_samples_count):
        coloring_samples = progeny_samples[i]
        assert coloring_samples in parents_samples
        mapping_from_progeny_to_paretns[str(i + 1)] = parents_samples.index(coloring_samples)
    mapping_from_progeny_to_paretns[str(coloring_samples_count + 1)] = -1
    return mapping_from_progeny_to_paretns


def get_vcf_header(samples):
    vcf_header = '##fileformat=VCFv4.2\n' \
                 '##GENOMAGIC_QUERY=generated with merge_vcf_pair.py in genomagic_qa\n' \
                 '##FORMAT=<ID=CO1,Number=3,Type=Integer,Description="Haplotypes Similarity Based Coloring (RGB)">\n' \
                 '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Overrides the length (in bases). Optional. ' \
                 'It can be present in some or all the variants on the line">\n' \
                 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in ' \
                 'this record">\n' \
                 '##INFO=<ID=HS,Number=0,Type=Flag,Description="Haplotypes similarity variant type">\n' \
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
    for s in samples:
        vcf_header = "{}\t{}".format(vcf_header,s)
    vcf_header = "{}\n".format(vcf_header)
    return vcf_header


def write_output_file(progeny_vcf_file, parents_vcf_file, out_vcf_file_name):
    parents_ite = vi.VcfIterator(parents_vcf_file)
    progeny_ite = vi.VcfIterator(progeny_vcf_file)
    progeny_samples = progeny_ite.get_headers()[vcf_info_columns_num:]
    parents_samples = parents_ite.get_headers()[vcf_info_columns_num:]
    assert progeny_ite.has_next()
    assert parents_ite.has_next()
    progeny_line_parts = progeny_ite.get_next_line_splitted()
    parents_line_parts = parents_ite.get_next_line_splitted()
    haps_list_of_progeny = progeny_line_parts[4]
    mapping_from_progeny_to_paretns = get_mapping_from_progeny_to_parents(haps_list_of_progeny, progeny_samples,
                                                                          parents_samples)
    output_vcf = open(out_vcf_file_name, 'w')
    output_vcf.write(get_vcf_header(parents_samples))
    while progeny_line_parts is not None and parents_line_parts is not None:
        [chr1, start1, end1] = get_position_from_vcf_splitted_line(progeny_line_parts)
        [chr2, start2, end2] = get_position_from_vcf_splitted_line(parents_line_parts)

        # We assume that the locations are the same for both files
        assert progeny_line_parts[4] == haps_list_of_progeny, "{} is not {}".format(progeny_line_parts[4],
                                                                                    haps_list_of_progeny)
        assert chr1 == chr2
        assert start1 < end2
        assert start2 < end1
        my_start = max(start1, start2)
        my_end = min(end1, end2)
        info_part = get_line_info(chr1, my_start, my_end, parents_line_parts[4])
        data_part = get_line_data(progeny_line_parts, mapping_from_progeny_to_paretns, parents_line_parts,
                                  my_end - my_start + 1)
        output_vcf.write("{}\t{}\n".format(info_part, data_part))
        if end1 == end2:
            progeny_line_parts = progeny_ite.get_next_line_splitted()
            parents_line_parts = parents_ite.get_next_line_splitted()
        elif end1 < end2:
            assert progeny_ite.has_next()
            progeny_line_parts = progeny_ite.get_next_line_splitted()
        else:
            assert end2 < end1
            assert parents_ite.has_next()
            parents_line_parts = parents_ite.get_next_line_splitted()
    assert not progeny_ite.has_next(), "progeny_ite.has_next()\n{}\n{}".format(progeny_line_parts, parents_line_parts)
    assert not parents_ite.has_next(), "parents_ite.has_next()\n{}\n{}".format(progeny_line_parts, parents_line_parts)
    output_vcf.close()


progeny_vcf_file = '/prodslow/testing/ariel/genomagic_qa/PSG-16/progeny_sim'
parents_vcf_file = '/prodslow/testing/ariel/genomagic_qa/PSG-16/parents_sim_fixed'
out_vcf_file_name = '/prodslow/testing/ariel/genomagic_qa/PSG-16/out'
write_output_file(progeny_vcf_file, parents_vcf_file, out_vcf_file_name)