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
