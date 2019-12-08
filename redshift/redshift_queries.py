import psycopg2
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
logging.getLogger().setLevel(logging.INFO)

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
    logging.info("executing {}".format(query))
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


def get_median_length_of_hap_similarity(host, data_version, min_score):
    temp_tables = get_hap_sim_of_hap_samples_sub_talble_string(host, data_version)
    query = 'WITH {}\nSELECT MEDIAN(len) OVER () as MEDIAN FROM hap_sim2 WHERE similarity_score >= {} limit 1;'.format(temp_tables,min_score)
    a = get_all_results(host, query)
    return int(a[0][0])


def get_average_length_of_hap_similarity(host, data_version, min_score):
    temp_tables = get_hap_sim_of_hap_samples_sub_talble_string(host, data_version)
    query = 'WITH {}\nselect AVG(len) from hap_sim2 where similarity_score >= {};'.format(temp_tables,min_score)
    a = get_all_results(host, query)
    return int(a[0][0])


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


def total_similarity(host, data_version):
    out_name = '{}/similarity_length.csv'.format(os.getcwd())
    similarity_table = get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    sql_query = "SELECT SUM(end_position-start_position) AS len , sample1, sample2 FROM {} where " \
                "identical_by_state='t' group by sample1, sample2;".format(similarity_table)
    rows = get_all_results(host, sql_query)
    df = pd.DataFrame.from_records(rows, columns=['len','sample1','sample2'])
    df.to_csv(out_name, index=False)
    return df['len'].mean()


def get_samples_table_by_analysis_method_query_string(host, data_version, analysis_methods):
    assert len(analysis_methods) > 0
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    filter_str = "analysis_method = \'{}\'".format(analysis_methods[0])
    for a in analysis_methods[1:]:
        filter_str = '{} OR analysis_method = \'{}\''.format(filter_str, a)
    query = 'select sample_id,analysis_method from {} where {}'.format(samples_table, filter_str)
    return query


def get_arg_wgs_full_haps_count(host, data_version):
    arg_wgs_query = get_samples_table_by_analysis_method_query_string(host, data_version, ['applied_reference_genome', 'whole_genome_sequencing'])
    haplotypes_samples_table = get_table_name(host, data_version, 'HAPLOTYPE_SAMPLES')
    haplotypes_info_table = get_table_name(host, data_version, 'HAPLOTYPES_INFO')
    samples_table_as = 'arg_wgs_samples AS ({})'.format(arg_wgs_query)
    haplotypes_samples_table_arg_wgs_as = 'samples AS (SELECT haplotype_idx, arg_wgs_samples.sample_id FROM {})'.format(
        get_inner_join_str(haplotypes_samples_table, 'arg_wgs_samples', 'sample_id', 'sample_id'))
    haps_info_as = 'haps AS (SELECT haplotype_idx, chromosome FROM {})'.format(haplotypes_info_table)
    full_query = 'WITH {}, {}, {} SELECT count(*),sample_id from {} where chromosome > 0 group by sample_id;'.format(
        samples_table_as, haplotypes_samples_table_arg_wgs_as, haps_info_as, get_inner_join_str('samples', 'haps', 'haplotype_idx', 'haplotype_idx'))
    #print(full_query)
    rows = get_all_results(host, full_query)
    out_name = '{}/haps_per_sample.csv'.format(os.getcwd())
    df = pd.DataFrame(rows, columns =['haps count', 'sample id'])
    df.to_csv(out_name, index=False)
    print("saved arg_wgs_full_haps_count to {}".format(out_name))
    #return get_all_results(host, full_query)


host = "rndlab-genomagic-redshift.cl6ox83ermwm.us-east-1.redshift.amazonaws.com"
data_version = "maize_benchmark_test_fix_mkrs_919_01"
s = get_arg_wgs_full_haps_count(host, data_version)
#df = pd.DataFrame(s, columns =['haps count', 'sample id'])
#print(df)




#def get_full_hap_samples_list():

 #   query = ''
#"""
#WITH arg_wgs_samples AS (select sample_id,analysis_method from maize_benchmark_test_fix_mkrs_919_01_samples_view where analysis_method = 'applied_reference_genome' OR analysis_method = 'whole_genome_sequencing'),
#samples AS (SELECT haplotype_idx, arg_wgs_samples.sample_id FROM maize_benchmark_test_fix_mkrs_919_01_haplotype_samples_view INNER JOIN arg_wgs_samples ON maize_benchmark_test_fix_mkrs_919_01_haplotype_samples_view.sample_id=arg_wgs_samples.sample_id),
#haps AS (SELECT haplotype_idx, chromosome FROM maize_benchmark_test_fix_mkrs_919_01_haplotypes_info_view)
#SELECT count(*),sample_id from samples INNER JOIN haps ON samples.haplotype_idx=haps.haplotype_idx where chromosome > 0 group by sample_id;"""

#def check_same_position():
#    'select count(haplotype_idx),length,chromosome,start_position,end_position from syn_maize_v4_3_gbs_haplotypes_info_view where chromosome=2 and start_position=end_position group by length,chromosome,start_position,end_position having count(haplotype_idx)'



#len_histogram(host, data_version)
#a = get_average_length_of_hap_similarity(host, data_version, 0)
#print(a)