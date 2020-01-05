import psycopg2
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
#logging.getLogger().setLevel(logging.INFO)

def create_connection_cursor(host):
    dbname='genomagic'
    user='nrgene'
    password=os.getenv("PGPASSWORD")
    assert password is not None, "env variable PGPASSWORD is not defined"
    port='5439'
    conn = psycopg2.connect("dbname={} user={} host={} port={} password={}".format(dbname,user,host,port,password))
    cur = conn.cursor()
    return cur


def get_all_results(host, query):
    assert query[-1] == ';'
    logging.info("executing sql query: {}".format(query))
    cur = create_connection_cursor(host)
    cur.execute(query)
    rows = cur.fetchall()
    cur.close()
    return rows


def get_table_name(host, data_version, table_type):
    query = 'SELECT table_name FROM {}_data_version WHERE table_type = \'{}\';'.format(data_version, table_type)
    rows = get_all_results(host, query)
    assert len(rows) == 1
    row = rows[0]
    assert len(row) == 1
    return row[0]


def get_sql_quey_as_data_frame(host, query, column_names):
    rows = get_all_results(host, query)
    df = pd.DataFrame(rows, columns=column_names)
    return df



def simple_grouping(host, table_name, groupby_field):
    query = 'SELECT COUNT(*),{} FROM {} GROUP BY {};'.format(groupby_field, table_name, groupby_field)
    return get_all_results(host, query)


def get_sample_types(host, data_version):
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    rows = simple_grouping(host, samples_table, 'analysis_method')
    df = pd.DataFrame(rows, columns =['count', 'analysis method'])
    return df


def printMarkersSummary(cur, data_version_name):
    print("summary of markers in {}:".format(data_version_name))
    cur.execute("""SELECT COUNT(*) FROM {}_markers;""".format(data_version_name))
    rows = cur.fetchall()
    alleles_count = int(rows[0][0])
    print("there are {} snp markers in the data version".format(alleles_count / 2))

    cur.execute("""SELECT count(CASE WHEN is_informative THEN 1 END) FROM {}_markers;""".format(data_version_name))
    rows = cur.fetchall()
    informative_count = int(rows[0][0])
    print("{} out of {} alleles are informative".format(informative_count, alleles_count))
    print("")


def printSamplesSummary(cur, data_version_name):
    print("summary of samples in {}:".format(data_version_name))
    analysis_methods = ["applied_reference_genome", "whole_genome_sequencing", "genotyping_by_sequencing", "snp_marker"]
    for a in analysis_methods:
        cur.execute(
            """SELECT COUNT(sample_idx),is_top_level FROM {}_samples WHERE analysis_method='{}' GROUP BY is_top_level;""".format(
                data_version_name, a))
        rows = cur.fetchall()
        n = len(rows)
        assert n <= 2
        if len(rows) > 0:
            x = [0, 0]
            for r in rows:
                x[int(r[1])] = int(r[0])
            print('{} {} samples, {} are top level'.format(sum(x), a, x[1]))
    print("")


def get_hap_info(cur, table_name, my_filter):
    cur.execute("""SELECT haplotype_idx,chromosome FROM {} WHERE {};""".format(table_name, my_filter))
    rows = cur.fetchall()
    print(rows)


def get_samples_type_info_as_string(host, data_version):
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    rows = simple_grouping(host, samples_table, 'analysis_method')
    s = ''
    for r in rows:
        assert len(r) == 2
        s = '{}\n{} samples count = {}'.format(s,r[1], r[0])
    return s


def get_hap_count_total_hap_markers_as_string(host, data_version):
    haplotypes_info_table = get_table_name(host, data_version, 'HAPLOTYPES_INFO')
    rows = get_all_results(host, 'SELECT COUNT(*) FROM {};'.format(haplotypes_info_table))
    hap_markers_count = int(rows[0][0])
    rows = get_all_results(host, 'SELECT COUNT(*) FROM {} WHERE chromosome=0;'.format(haplotypes_info_table))
    unmapped_hap_markers_count = int(rows[0][0])
    return [hap_markers_count, unmapped_hap_markers_count]
    #'in table {} there are total of {}M haplotype markers {}M of them are unmapped'.format(haplotypes_info_table, hap_markers_count, unmapped_hap_markers_count)


def get_hap_samples_total_as_string(host, data_version):
    haplotypes_samples_table = get_table_name(host, data_version, 'HAPLOTYPE_SAMPLES')
    sub_table_samples = 'samples AS (SELECT haplotype_idx, sample_id FROM {})'.format(haplotypes_samples_table)
    haplotypes_info_table = get_table_name(host, data_version, 'HAPLOTYPES_INFO')
    sub_table_haps = 'haps AS (SELECT haplotype_idx, chromosome FROM {})'.format(haplotypes_info_table)
    inner_join_sub_table = 'samples INNER JOIN haps ON samples.haplotype_idx=haps.haplotype_idx'
    hapsXsamples_query = 'WITH {}, {} SELECT COUNT(*) FROM {}'.format(sub_table_samples, sub_table_haps,
                                                                      inner_join_sub_table)
    rows = get_all_results(host, '{};'.format(hapsXsamples_query))
    total_haps_samples = int(rows[0][0])
    rows = get_all_results(host, '{} WHERE chromosome=0;'.format(hapsXsamples_query))
    unmapped_haps_samples = int(rows[0][0])
    return [total_haps_samples, unmapped_haps_samples]
    #print('There are {}M hapsXsamples {}M of them are unmapped'.format(total_haps_samples, unmapped_haps_samples))


def hist_count(host, data_version):
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    score_query = 'SELECT COUNT(*), ROUND(similarity_score,1) as sc from {} GROUP BY sc;'.format(similarity_table)

    rows = get_all_results(host, score_query)
    assert len(rows) == 11
    df = pd.DataFrame.from_records(rows, columns=['Count', 'Score']).sort_values(by=['Score'])
    height = df['Count']
    bars = df['Score']
    y_pos = np.arange(len(bars))
    plt.figure()
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.show()


def hist_count2(host, data_version):
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    score_len_query = 'SELECT AVG(end_position-start_position), ROUND(similarity_score,2) as sc from {} GROUP BY sc;'.format(similarity_table)
    rows = get_all_results(host, score_len_query)
    df = pd.DataFrame.from_records(rows, columns =['mean len', 'score']).sort_values(by=['score'])
    plt.figure()
    plt.scatter(df['score'], df['mean len'])
    plt.title('score vs len')
    plt.xlabel('score')
    plt.ylabel('len')
    plt.show()


def get_inner_join_str(table1, table2, field1, field2):
    return '{} INNER JOIN {} ON {}.{}={}.{}'.format(table1, table2, table1, field1, table2, field2)


def get_hap_sim_of_hap_samples_sub_talble_string(host, data_version):
    """
    This is a function that returns the haplotype similarity , filtered on wgs/arg/gbs samples - so we exclude snp similarity
    :return: the string to create the filtered haplotype similarity in the sql
    """
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    samples_with_haps = 'analysis_method=\'applied_reference_genome\' OR analysis_method=\'whole_genome_sequencing\' OR analysis_method=\'genotyping_by_sequencing\''
    arg_wgs_sub_table = 'wgs_samples AS (SELECT sample_id FROM {} WHERE  {})'.format(samples_table, samples_with_haps)
    inner_join_on_sample1 = get_inner_join_str ('wgs_samples', similarity_table, 'sample_id', 'sample1')
    hap_sim_sample1_wgs = 'hap_sim1 AS (SELECT * from {})'.format(inner_join_on_sample1)
    inner_join_on_sample2 = get_inner_join_str ('wgs_samples', 'hap_sim1', 'sample_id', 'sample2')
    hap_sim_sample2_wgs = 'hap_sim2 AS (SELECT similarity_score, end_position-start_position as len FROM {})'.format(inner_join_on_sample2)
    return '{},\n{},\n{}'.format(arg_wgs_sub_table, hap_sim_sample1_wgs, hap_sim_sample2_wgs)





def len_histogram(host, data_version):
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    sim_query = 'SELECT COUNT(*), ROUND(log(end_position-start_position),1) as sc from {} GROUP BY sc ORDER BY sc;'.format(similarity_table)
    rows = get_all_results(host, sim_query)
    df = pd.DataFrame.from_records(rows, columns=['Count', 'Length']).sort_values(by=['Length'])
    height = df['Count']
    bars = df['Length']
    y_pos = np.arange(len(bars))
    plt.figure()
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.show()





def get_samples_table_by_analysis_method_query_string(host, data_version, analysis_methods):
    assert len(analysis_methods) > 0
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    filter_str = "analysis_method = \'{}\'".format(analysis_methods[0])
    for a in analysis_methods[1:]:
        filter_str = '{} OR analysis_method = \'{}\''.format(filter_str, a)
    query = 'select sample_id,analysis_method from {} where {}'.format(samples_table, filter_str)
    return query

def get_mapped_haplotype_samples_table_only_arg_wgs_string(host, data_version, temp_table_name1, temp_table_name2, temp_table_name3):
    table_1_content = get_samples_table_by_analysis_method_query_string(host, data_version, ['applied_reference_genome',
                                                                                           'whole_genome_sequencing'])
    haplotypes_samples_table = get_table_name(host, data_version, 'HAPLOTYPE_SAMPLES')
    haplotypes_info_table = get_table_name(host, data_version, 'HAPLOTYPES_INFO')
    inner_join_str = get_inner_join_str(haplotypes_samples_table, temp_table_name1, 'sample_id', 'sample_id')
    table_2_content = 'SELECT haplotype_idx, {}.sample_id, analysis_method FROM {}'.format(temp_table_name1, inner_join_str)
    inner_join_str2 = get_inner_join_str(haplotypes_info_table, temp_table_name2, 'haplotype_idx', 'haplotype_idx')
    table_3_content = 'SELECT sample_id,analysis_method,{}.haplotype_idx from {} WHERE chromosome > 0'.format(temp_table_name2, inner_join_str2)
    temp_tables_string = '{} as ({}), {} as ({}), {} as ({})'.format(temp_table_name1, table_1_content, temp_table_name2,
                                                                  table_2_content, temp_table_name3, table_3_content)
    return temp_tables_string


def write_samples_haps_count_to_file(host, data_version):
    mapped_table_name = 'mapped_hap_samples'
    haps_samples_mapped_arg_wgs = get_mapped_haplotype_samples_table_only_arg_wgs_string(host, data_version,'temp_table1', 'temp_table2', mapped_table_name)
    full_query = 'WITH {} SELECT sample_id, analysis_method, count(haplotype_idx) FROM {} GROUP BY sample_id, analysis_method;'.format(haps_samples_mapped_arg_wgs, mapped_table_name)
    rows = get_all_results(host, full_query)
    out_name = '{}/{}_haps_per_sample.csv'.format(os.getcwd(), data_version)
    df = pd.DataFrame(rows, columns =['sample id', 'analysis_method', 'haps count'])
    df.to_csv(out_name, index=False)
    print("saved arg_wgs_full_haps_count to {}".format(out_name))


def flat_freq_dataframe_to_matrix(df):
    max_freq = df['hap freq'].max()
    i = 1
    temp_df1 = df.loc[df['hap freq'] == i, ['sample id', 'count']].rename(columns={'count': 'freq {}'.format(i)})
    for i in range(2, max_freq + 1):
        temp_df2 = df.loc[df['hap freq'] == i, ['sample id', 'count']].rename(columns={'count': 'freq {}'.format(i)})
        temp_df1 = pd.merge(temp_df1, temp_df2, on='sample id', how='outer')
    return temp_df1.fillna(0)


def write_samples_haps_freq_to_file(host, data_version):
    mapped_table_name = 'mapped_hap_samples'
    haps_freq_table_name = 'haps_freq'
    haps_samples_mapped_arg_wgs = get_mapped_haplotype_samples_table_only_arg_wgs_string(host, data_version,
                                                                                         'temp_table1', 'temp_table2',
                                                                                         mapped_table_name)
    haps_freq = '{} as (select count(*) as freq,haplotype_idx from {} group by haplotype_idx)'.format(haps_freq_table_name, mapped_table_name)
    inner_join_str = get_inner_join_str(haps_freq_table_name, mapped_table_name, 'haplotype_idx', 'haplotype_idx')
    full_query = 'WITH {}, {} SELECT sample_id,freq,count(haps_freq.haplotype_idx) FROM {} GROUP BY sample_id,freq ORDER BY sample_id,freq;'.format(
        haps_samples_mapped_arg_wgs, haps_freq, inner_join_str)
    rows = get_all_results(host, full_query)
    out_name = '{}/{}_haps_freq_per_sample.csv'.format(os.getcwd(), data_version)
    df = pd.DataFrame(rows, columns=['sample id', 'hap freq', 'count'])
    mat_df = flat_freq_dataframe_to_matrix(df)
    mat_df.to_csv(out_name, index=False)
    print("saved arg_wgs_full_haps_count to {}".format(out_name))
    return df


def get_haplotype_similarity_with_analysis_types(host, data_version, table_name):
    temp_var1 = 'temp_var1'
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    fields1 = 'sample1, analysis_method as sample1_type,sample2,end_position-start_position as len, similarity_score'
    inner_join1 = get_inner_join_str(similarity_table, samples_table, 'sample1', 'sample_id')
    table1 = '{} AS (SELECT {} FROM {})'.format(temp_var1, fields1, inner_join1)
    fields2 = 'sample1, sample1_type,sample2, analysis_method as sample2_type, len, similarity_score'
    inner_join2 = get_inner_join_str(temp_var1, samples_table, 'sample2', 'sample_id')
    table2 = '{} AS (SELECT {} FROM {})'.format(table_name, fields2, inner_join2)
    return '{}, {}'.format(table1, table2)


def get_specific_sample_types_string(pairs_list, field1, field2):
    my_str = '({}=\'{}\' AND {}=\'{}\')'.format(field1, pairs_list[0][0], field2,pairs_list[0][1])
    for p in pairs_list[1:]:
        my_str = '{} OR {}'.format(my_str, '({}=\'{}\' AND {}=\'{}\')'.format(field1, p[0], field2, p[1]))
    return my_str


def get_wgs_haplotype_similarity(host, data_version, pairs_list, table_name):
    pairs_str = get_specific_sample_types_string(pairs_list, 'sample1_type', 'sample2_type')
    hap_sim_types_name = 'hap_sim_types_table'
    hap_sim_types_str = get_haplotype_similarity_with_analysis_types(host, data_version, hap_sim_types_name)
    query = '{}, {} as (SELECT len,similarity_score FROM {} WHERE {})'.format(hap_sim_types_str, table_name, hap_sim_types_name, pairs_str)
    return query



def get_median_length_of_hap_similarity(host, data_version, pairs_list, min_score):
    temp_table_name = 'score_and_len_from_hap_sim_table'
    temp_tables = get_wgs_haplotype_similarity(host, data_version, pairs_list, temp_table_name)
    query = 'WITH {} SELECT MEDIAN(len) OVER () FROM {} WHERE similarity_score >= {} LIMIT 1;'.format(temp_tables, temp_table_name, min_score)
    a = get_all_results(host, query)
    assert a is not None
    if len(a) == 0:
        return None
    assert a[0] is not None
    assert a[0][0] is not None
    return int(a[0][0])


def get_average_length_of_hap_similarity(host, data_version, pairs_list, min_score):
    temp_table_name = 'score_and_len_from_hap_sim_table'
    temp_tables = get_wgs_haplotype_similarity(host, data_version, pairs_list, temp_table_name)
    query = 'WITH {} SELECT AVG(len) FROM {} where similarity_score >= {};'.format(temp_tables, temp_table_name, min_score)
    a = get_all_results(host, query)
    assert a is not None
    assert len(a) > 0
    assert a[0] is not None
    if a[0][0] is None:
        return None
    else:
        return int(a[0][0])


def compute_sim_len_histogram(host, data_version, pairs_list, min_score):
    temp_table_name = 'score_and_len_from_hap_sim_table'
    temp_tables = get_wgs_haplotype_similarity(host, data_version, pairs_list, temp_table_name)
    query = 'WITH {} SELECT COUNT(*), ROUND(log(len),1) as sc FROM {} where similarity_score >= {} GROUP BY sc ORDER BY sc;'.format(temp_tables, temp_table_name,
                                                                                   min_score)

    rows = get_all_results(host, query)
    df = pd.DataFrame.from_records(rows, columns=['Count', 'Length']).sort_values(by=['Length'])
    height = df['Count']
    bars = df['Length']
    y_pos = np.arange(len(bars))
    plt.figure()
    plt.bar(y_pos, height)
    plt.xticks(y_pos, bars)
    plt.show()


def write_all_pairwise_similarities(host, data_version, pairs_list, threshold, out_name):
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    temp_table1 = 'temp_var1 AS (SELECT sample1, analysis_method as sample1_type,sample2,' \
                  'end_position-start_position as len, similarity_score FROM {} INNER JOIN {} ON ' \
                  '{}.sample1={}.sample_id)'.format(similarity_table, samples_table, similarity_table, samples_table)
    temp_table2 = 'temp_var2 AS (SELECT sample1, sample1_type,sample2, analysis_method as sample2_type, len, ' \
                  'similarity_score FROM temp_var1 INNER JOIN {} ON ' \
                  'temp_var1.sample2={}.sample_id)'.format(samples_table, samples_table)
    pairs_str = get_specific_sample_types_string(pairs_list, 'sample1_type', 'sample2_type')
    temp_table3  = 'temp_var3 as (SELECT  len , similarity_score,sample1, sample2 FROM temp_var2 where {})'.format(pairs_str)
    query = 'WITH {},{},{} SELECT sample1, sample2, SUM(len) FROM temp_var3 ' \
            'WHERE similarity_score>={} GROUP BY sample1, sample2;'.format(temp_table1, temp_table2, temp_table3, threshold)
    #out_name = '{}/similarity_length_per_comparison.csv'.format(os.getcwd())
    rows = get_all_results(host, query)
    df = pd.DataFrame.from_records(rows, columns=['sample1','sample2', 'total similarity length'])
    print('average similarity length of comparisons is {}'.format(df['total similarity length'].mean()))
    df.to_csv(out_name, index=False)


def get_genome_size(host, data_version):
    pivot_name = get_table_name(host, data_version, 'PIVOT_SAMPLE')
    chr_len_table = get_table_name(host, data_version, 'CHROMOSOME_LENGTH')
    query = 'SELECT SUM(chromosome_length) FROM {} where sample=\'{}\';'.format(chr_len_table, pivot_name)
    rows = get_all_results(host, query)
    genome_size = int(rows[0][0])
    return genome_size



#data_version = 'public_soy_v2_03'
#host='rndlab-genomagic-redshift.cl6ox83ermwm.us-east-1.redshift.amazonaws.com'
#get_genome_size(host, data_version)



#sample_pair_types = []
#sample_pair_types.append(['whole_genome_sequencing', 'whole_genome_sequencing'])
#sample_pair_types.append(['applied_reference_genome', 'whole_genome_sequencing'])
#

#my_query = 'WITH wgs_samples as (select sample_id from public_soy_v2_03_samples_view where analysis_method=\'applied_reference_genome\' or analysis_method=\'whole_genome_sequencing\'), sim_1 as (select sample1,sample2,chromosome_id,start_position,end_position,similarity_score,identical_by_state from public_soy_v2_03_haplotypes_similarity_view inner join wgs_samples on public_soy_v2_03_haplotypes_similarity_view.sample1=wgs_samples.sample_id) , sim_2 as (select * from sim_1 inner join wgs_samples on sim_1.sample2=wgs_samples.sample_id) select sample1,sample2,chromosome_id,start_position,end_position,similarity_score,identical_by_state from sim_2 where identical_by_state=\'t\' and chromosome_id=1 and start_position<=10000000 and end_position>=10000000;'




#pairs_list = []
#pairs_list.append(['whole_genome_sequencing', 'whole_genome_sequencing'])
#pairs_list.append(['applied_reference_genome', 'whole_genome_sequencing'])
#pairs_list.append(['whole_genome_sequencing', 'applied_reference_genome'])
#pairs_list.append(['applied_reference_genome', 'applied_reference_genome'])
#
#temp_table_name = 'score_and_len_from_hap_sim_table'
#

#sim_table = get_haplotype_similarity_with_analysis_types(host, data_version, temp_table_name)
#print(sim_table)
#pairs_str = get_specific_sample_types_string(pairs_list, 'sample1_type', 'sample2_type')
#print(pairs_str)
#print("WITH {}, arg_wgs_sim as (select * from {} where {}) select * from arg_wgs_sim limit 10;".format(sim_table, temp_table_name, pairs_str))
#print(pairs_str)

#data_version='maize_benchmark_test_fix_mkrs_919_01'
#write_all_pairwise_similarities(host, data_version, pairs_list, 0.99)
#s = get_average_length_of_hap_similarity(host, data_version, pairs_list, 0)
#print (s)

