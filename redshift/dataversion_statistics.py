from . import redshift_basic as rs
import pandas as pd


def get_sample_types(host, data_version):
    samples_table = rs.get_table_name(host, data_version, 'SAMPLES')
    rows = rs.simple_grouping(host, samples_table, 'analysis_method')
    return rows


