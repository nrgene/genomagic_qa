import os
import sys
genomagic_abs_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append('{}/vcf_utils'.format(genomagic_abs_path))
import vcf_iterator as vi


def return_tsv_line(sample_name, marker_id, alleles, flag):
    assert len(alleles) == 2, '{} {} {}'.format(alleles,marker_id, flag)
    if flag == 0:
        return ''
    status1, status2 = get_status_from_flag(flag)
    line1 = '{}\t{}\t{}\t{}\n'.format(sample_name, marker_id, status1, alleles[0])
    line2 = '{}\t{}\t{}\t{}\n'.format(sample_name, marker_id, status2, alleles[1])
    return '{}{}'.format(line1, line2)


def get_status_from_flag(flag):
    if flag == 1:
        return 1, 0
    elif flag == 2:
        return 0, 1
    assert flag == 3, flag
    return 1, 1


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


def handle_vcf_line(splitted_vcf_line, vcf_columns_headers, marker_data_from_vcf_line, genotype_to_flag, columns_setoff):
    assert len(splitted_vcf_line) == len(vcf_columns_headers), splitted_vcf_line
    [marker_id, alleles] = marker_data_from_vcf_line(splitted_vcf_line)
    my_str = ''
    columns_count = len(vcf_columns_headers)
    for i in range(columns_setoff, columns_count):
        sample_name = vcf_columns_headers[i]
        current_genotype = splitted_vcf_line[i]
        tsv_line = get_tsv_line_from_genotype(alleles, current_genotype, genotype_to_flag, marker_id, sample_name)
        my_str = '{}{}'.format(my_str, tsv_line)
    return my_str


def get_tsv_line_from_genotype(alleles, current_genotype, genotype_to_flag, marker_id, sample_name):
    assert current_genotype in genotype_to_flag
    flag = genotype_to_flag[current_genotype]
    tsv_line = return_tsv_line(sample_name, marker_id, alleles, flag)
    return tsv_line


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
    f = open(output_file, "w")
    f.write('sample_id\tmarker_id\tstatus\tallele_seq\n')
    while vct_ite.has_next():
        curr_line = vct_ite.get_next_line_splitted()
        tsv_lines = handle_vcf_line(curr_line, headers, marker_data_from_vcf_line, genotype_to_flag, columns_setoff)
        output_file.write(tsv_lines)
    f.close()
