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

def get_all_results(host, query):
    cur = createConnectionCursor(host)
    cur.execute(query)
    rows = cur.fetchall()
    cur.close()
    return rows


def get_table_name(host, data_version, table_name):
    query = 'select table_name from {}_data_version where table_type = \'{}\';'.format(data_version, table_name)
    cur = createConnectionCursor(host)
    cur.execute(query)
    row = cur.fetchone()
    cur.close()
    return row




