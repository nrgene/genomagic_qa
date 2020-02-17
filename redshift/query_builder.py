def get_inner_join_str(table1, table2, field1, field2):
    return '{} INNER JOIN {} ON {}.{}={}.{}'.format(table1, table2, table1, field1, table2, field2)


def get_specific_sample_types_string(pairs_list, field1, field2):
    my_str = '({}=\'{}\' AND {}=\'{}\')'.format(field1, pairs_list[0][0], field2,pairs_list[0][1])
    for p in pairs_list[1:]:
        my_str = '{} OR {}'.format(my_str, '({}=\'{}\' AND {}=\'{}\')'.format(field1, p[0], field2, p[1]))
    return my_str


def get_all_pairwise_similarities_query(pairs_list, threshold, similarity_table_name, samples_table_name):
    temp_table1 = 'temp_var1 AS (SELECT sample1, analysis_method as sample1_type,sample2,' \
                  'end_position-start_position as len, similarity_score FROM {} INNER JOIN {} ON ' \
                  '{}.sample1={}.sample_id)'.format(similarity_table_name, samples_table_name, similarity_table_name,
                                                    similarity_table_name)
    temp_table2 = 'temp_var2 AS (SELECT sample1, sample1_type,sample2, analysis_method as sample2_type, len, ' \
                  'similarity_score FROM temp_var1 INNER JOIN {} ON ' \
                  'temp_var1.sample2={}.sample_id)'.format(samples_table_name, samples_table_name)
    pairs_str = get_specific_sample_types_string(pairs_list, 'sample1_type', 'sample2_type')
    temp_table3 = 'temp_var3 as (SELECT  len , similarity_score,sample1, sample2 FROM temp_var2 where {})'.format(
        pairs_str)
    query = 'WITH {},{},{} SELECT sample1, sample2, SUM(len) FROM temp_var3 ' \
            'WHERE similarity_score>={} GROUP BY sample1, sample2;'.format(temp_table1, temp_table2, temp_table3,
                                                                           threshold)
    return query


def get_haplotype_similarity_with_analysis_types( table_name, similarity_table, samples_table):
    temp_var1 = 'temp_var1'
    fields1 = 'sample1, analysis_method as sample1_type,sample2,end_position-start_position as len, similarity_score'
    inner_join1 = get_inner_join_str(similarity_table, samples_table, 'sample1', 'sample_id')
    table1 = '{} AS (SELECT {} FROM {})'.format(temp_var1, fields1, inner_join1)
    fields2 = 'sample1, sample1_type,sample2, analysis_method as sample2_type, len, similarity_score'
    inner_join2 = get_inner_join_str(temp_var1, samples_table, 'sample2', 'sample_id')
    table2 = '{} AS (SELECT {} FROM {})'.format(table_name, fields2, inner_join2)
    return '{}, {}'.format(table1, table2)


def get_length_score_by_sample_types(pairs_list, table_name, similarity_table, samples_table):
    pairs_str = get_specific_sample_types_string(pairs_list, 'sample1_type', 'sample2_type')
    hap_sim_types_name = 'hap_sim_types_table'
    hap_sim_types_str = get_haplotype_similarity_with_analysis_types( hap_sim_types_name, similarity_table, samples_table)
    query = '{}, {} AS (SELECT len,similarity_score FROM {} WHERE {})'.format(hap_sim_types_str, table_name, hap_sim_types_name, pairs_str)
    return query


def get_haplotype_frequency(hap_samples_table, hap_info, samples_table):
    sample_types = ['whole_genome_sequencing', 'whole_genome_sequencing']
    wgs_samples_sub_query = get_samples_filtered_by_analysis_types(samples_table, sample_types)
    # step1: sub_table1 is the arg/wgs samples
    sub_table1 = 'temp_var1 AS ({})'.format(wgs_samples_sub_query)
    # step2: sub_table2 is hap samples table only for the arg/wgs samples
    sub_table2 = 'temp_var2 AS (SELECT haplotype_idx,temp_var1.sample_id FROM {})'.format(
        get_inner_join_str('temp_var1', hap_samples_table, 'sample_id', 'sample_id'))
    # step3: sub_table3 are the mapped haplotypes
    sub_table3 = 'temp_var3 AS (SELECT haplotype_idx FROM {} WHERE chromosome > 0)'.format(hap_info)
    # step4: sub_table4 are the mapped haplotypes and the arg/wgs samples
    sub_table4 = 'temp_var4 AS (SELECT sample_id,temp_var2.haplotype_idx FROM {})'.format(
        get_inner_join_str('temp_var2', 'temp_var3', 'haplotype_idx', 'haplotype_idx'))
    # step5: sub_table4 are the mapped haplotypes and the arg/wgs samples
    sub_table5 = 'temp_var5 AS (SELECT haplotype_idx,count(sample_id) as freq FROM temp_var4 GROUP BY haplotype_idx)'
    sub_table6 = 'temp_var6 AS (SELECT freq,count(haplotype_idx) FROM temp_var5 GROUP BY freq)'
    merged_sub_tables = ', '.join([sub_table1, sub_table2, sub_table3, sub_table4, sub_table5, sub_table6])
    query ='WITH {} SELECT freq,count FROM temp_var6 ORDER BY freq;'.format(merged_sub_tables)
    return query



def get_samples_filtered_by_analysis_types(samples_table, analysis_types_list):
    expressions_list = ['analysis_method=\'{}\''.format(s) for s in analysis_types_list]
    filter_string = ' OR '.join(expressions_list)
    return 'SELECT * FROM {} WHERE {}'.format(samples_table, filter_string)