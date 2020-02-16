import sys
import numpy as np
genomagic_qa_repo_path = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(genomagic_qa_repo_path)
import redshift.basic_queries as rs

def get_similarities_in_position(host, data_version, chromosome, position, threshold):
    query = 'select * from public_soy_v1_16_haplotypes_similarity_view where sample1 = \'s_4j105_3_4\' ' \
            'and identical_by_state = \'t\' and chromosome_id = 1 and start_position <= 14000000 and end_position >= 14000000;'


