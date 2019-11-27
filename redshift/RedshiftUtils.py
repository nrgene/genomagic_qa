import psycopg2
import os

def createConnectionCursor(host):
    dbname='genomagic'
    user='nrgene'
    password=os.getenv("PGPASSWORD")
    assert password is not None, "env variable PGPASSWORD is not defined"
    port='5439'
    conn = psycopg2.connect("dbname={} user={} host={} port={} password={}".format(dbname,user,host,port,password))
    cur = conn.cursor()
    return cur


def printMarkersSummary(cur,data_version_name):
    print("summary of markers in {}:".format(data_version_name))
    cur.execute("""SELECT COUNT(*) FROM {}_markers;""".format(data_version_name))
    rows = cur.fetchall()
    alleles_count=int(rows[0][0])
    print("there are {} snp markers in the data version".format(alleles_count/2))
    
    cur.execute("""SELECT count(CASE WHEN is_informative THEN 1 END) FROM {}_markers;""".format(data_version_name))
    rows = cur.fetchall()
    informative_count=int(rows[0][0])
    print("{} out of {} alleles are informative".format(informative_count,alleles_count))
    print("")

def printSamplesSummary(cur,data_version_name):
    print("summary of samples in {}:".format(data_version_name))
    analysis_methods=["applied_reference_genome", "whole_genome_sequencing", "genotyping_by_sequencing", "snp_marker"]
    for a in analysis_methods:
        cur.execute("""SELECT COUNT(sample_idx),is_top_level FROM {}_samples WHERE analysis_method='{}' GROUP BY is_top_level;""".format(data_version_name,a))
        rows = cur.fetchall()
        n=len(rows)
        assert n<=2
        if len(rows)>0:
            x=[0,0]
            for r in rows:
                x[int(r[1])]=int(r[0])
            print('{} {} samples, {} are top level'.format(sum(x),a,x[1]))
    print ("")


def get_hap_info(cur,table_name, my_filter):
    cur.execute("""SELECT haplotype_idx,chromosome FROM {} WHERE {};""".format(table_name, my_filter))
    rows = cur.fetchall()
    print(rows)
    #alleles_count=int(rows[0][0])
    #print "there are {} snp markers in the data version".format(alleles_count/2)

