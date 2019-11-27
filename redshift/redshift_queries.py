import psycopg2
import os
import pandas as pd

def createConnectionCursor(host):
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
    cur = createConnectionCursor(host)
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


def get_samples_type_info(host, data_version):
    samples_table = get_table_name(host, data_version, 'SAMPLES')
    rows = simple_grouping(host, samples_table, 'analysis_method')
    s = ''
    for r in rows:
        assert len(r) == 2
        print('{} samples count = {}'.format(r[1], r[0]))

def
haplotypes_info_table = rs.get_table_name(host, data_version, 'HAPLOTYPES_INFO')
rows = rs.get_all_results(host, 'SELECT COUNT(*) FROM {};'.format(haplotypes_info_table))
hap_markers_count = int(rows[0][0]/1000000)
rows = rs.get_all_results(host, 'SELECT COUNT(*) FROM {} WHERE chromosome=0;'.format(haplotypes_info_table))
unmapped_hap_markers_count = int(rows[0][0]/1000000)
print('in table {} there are total of {}M haplotype markers {}M of them are unmapped'.format(haplotypes_info_table, hap_markers_count, unmapped_hap_markers_count))











