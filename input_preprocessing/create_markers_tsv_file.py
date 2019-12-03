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



my_file = '/prodslow/testing/ariel/genomagic_qa/PSG-16/SoyNAM_parents+progeny_4312_SNP_Wm82.a2.filtered'
vct_ite = VcfIterator(my_file)
assert vct_ite.has_next()
line = vct_ite.get_next_line()
print(line)


