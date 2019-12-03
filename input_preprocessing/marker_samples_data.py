import os
import sys
genomagic_abs_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append('{}/vcf_utils'.format(genomagic_abs_path))
import vcf_iterator as vi

def return_tsv_line(sample_name, marker_id, alleles, flag):
    assert len(alleles) == 2, '{} {} {}'.format(alleles,marker_id, flag)
    status1 = 0
    status2 = 0
    if flag == 0:
        return ''
    elif flag == 1:
        status1 = 1
    elif flag == 2:
        status2 = 1
    else:
        assert flag == 3, flag
        status1 = 1
        status2 = 1
    line1 = '{}\t{}\t{}\t{}\n'.format(sample_name, marker_id, alleles[0], status1)
    line2 = '{}\t{}\t{}\t{}\n'.format(sample_name, marker_id, alleles[1], status2)
    return '{}{}'.format(line1, line2)


def vcf_line_to_tsv_lines(vcf_splitted_line, headers, alleles, marker_id, genotype_to_flag, columns_setoff):
    my_str = ''
    columns_count = len(headers)
    for i in range(columns_setoff, columns_count):
        sample_name = headers[i]
        assert vcf_splitted_line[i] in genotype_to_flag
        flag = genotype_to_flag[vcf_splitted_line[i]]
        tsv_line = return_tsv_line(sample_name, marker_id, alleles, flag)
        my_str = '{}{}'.format(my_str, tsv_line)
    return my_str


def write_tsv_file(input_file, output_file, genotype_to_flag, marker_data_from_vcf_line, columns_setoff):
    """The main function to call from outside that creates the tsv file
    parameters:
    input_file - the input vcf file

    output_file - the tsv outfile that will be created

    genotype_to_flag - dict that gets allele and returns flag (0=none,1=first allele,2=second allele,3=both alleles)
    example : my_dict = {'-': 0, 'A': 1, 'B': 2, 'H':3}

    marker_data_from_vcf_line - function that gets the splitted vcf line and returns marker_id and alleles

    columns_setoff - how many columns are info columns
    """
    vct_ite = vi.VcfIterator(input_file)
    assert vct_ite.has_next()
    headers = vct_ite.get_headers()
    headers_num = len(headers)
    f = open(output_file, "w")
    f.write('sample_id\tmarker_id\tallele_seq\tstatus\n')
    while vct_ite.has_next():
        curr_line = vct_ite.get_next_line_splitted()
        assert len(curr_line) == headers_num, curr_line
        [marker_id, alleles] = marker_data_from_vcf_line(curr_line)
        tsv_lines = vcf_line_to_tsv_lines(curr_line, headers, alleles, marker_id, genotype_to_flag, columns_setoff)
        f.write(tsv_lines)
    f.close()
