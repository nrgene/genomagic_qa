import sys
genomagic_path = '/home/ariel/clients/genomagic_qa/'
sys.path.append(genomagic_path)
import redshift.redshift_queries as rq
import api_server.requests as api
import analysis.markers_comparison as mc
import analysis.similarities_comparison as sc
"""
In this script we are extracting snp data from wgs samples, and then generate similarities based on these snps
These similarities are compared to haplotype based similarities from redshift
The assumption is that this is a way to optimize both hap similarity parameters and markers similarity parameteres
"""


def test_sim_by_data_version(host, data_version, snp_sim, samples_list):
    similarity_table = rq.get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    query = 'select sample1, sample2, chromosome_id, start_position, end_position from {} where chromosome_id={} and similarity_score>={} order by start_position;'.format(similarity_table, chromosome, similarity_threshold)
    hap_similarities = rq.get_sql_query_as_data_frame(host, query, ['s1','s2','chr','start','end'])
    print("hap_similarities shape = {}".format(hap_similarities.shape))
    df = sc.get_similarities_comparison_by_samples_pair(snp_sim, hap_similarities, samples_list)
    unique_to_snp = df['unique in 1'].sum()
    unique_to_hap = df['unique in 2'].sum()
    shared_len = df['both'].sum()
    return [unique_to_snp, unique_to_hap, shared_len]


host='rndlab-genomagic-redshift.cl6ox83ermwm.us-east-1.redshift.amazonaws.com'
api_server='api-dev.nrgene.local'
samples = 'lee__ver100,pi483463__ver100,williams82__ver100,fc33243_anderson__ver100,pi180501_strain_no__18__ver100,pi240664_bilomi_no__3__ver100,pi266806c_no_4__ver100,pi322692_max_c_p1159a8__ver100,pi360957_karafuto_no__1__ver100,pi361087_medias_23__ver100,pi361093_novosadska_br__1__ver100,pi437265d_dobruzanca_d__ver100,pi438335_sao_196_c__ver100,pi468908__ver100,pi497964a_i_c__9461__ver100,pi507180_rikuu_21__ver100,pi518668_tn_4_86__ver100,pi548193_t201__ver100,pi548360_korean__ver100,pi548364_macoupin__ver100,pi548447_cherokee__ver100,pi548452_dixie__ver100,pi548490_tanner__ver100,pi548520_preston__ver100,pi548561_hodgson__ver100,pi548619_sparks__ver100,pi549041a_zyd_2709__ver100,pi559932_manokin__ver100,pi567426_bai_huang_dou__ver100,pi567558_liu_shi_ri_jin_huang_da_dou__ver100,pi567604a_xin_huang_dou__ver100,pi567788_bienville__ver100,pi578309_i_64__ver100,pi592523_glacier__ver100'
max_hap_num = 1000000
similarity_threshold=0.85
chromosome=2
min_p = 1
win_len = 20
trim_len = 0
min_samples_presence = 2
max_major_allele_freq = 0.8
samples_list = samples.split(',')

data_version = 'public_soy_v1_14' # should not afect the similarity
snp_markers = api.get_snp_markers_as_dataframe(api_server, data_version, samples, chromosome, max_hap_num)
snp_sim = mc.get_all_pairwise_similarities_snp(snp_markers, samples_list, chromosome, min_p, win_len, trim_len, max_major_allele_freq, min_samples_presence)


[a1,a2,a3] = test_sim_by_data_version(host, 'public_soy_v1_14', snp_sim, samples_list)
print("1.0 jacard = {}".format(a3 / (a3 +  a1 + a2)))
print("missing similarities fraction={}".format(a1/(a1+a3)))
print("extra similarities fraction={}".format(a2/(a2+a3)))

[b1,b2,b3] = test_sim_by_data_version(host, 'public_soy_test_remove_freq', snp_sim, samples_list)
print("0.9 jacard = {}".format(b3 / (b3 +  b1 + b2)))
print("missing similarities fraction={}".format(b1/(b1+b3)))
print("extra similarities fraction={}".format(b2/(b2+b3)))

[c1,c2,c3] = test_sim_by_data_version(host, 'public_soy_test_remove_freq_3', snp_sim, samples_list)
print("0.8 jacard = {}".format(c3 / (c3 +  c1 + c2)))
print("missing similarities fraction={}".format(c1/(c1+c3)))
print("extra similarities fraction={}".format(c2/(c2+c3)))

[d1,d2,d3] = test_sim_by_data_version(host, 'public_soy_test_remove_freq_2', snp_sim, samples_list)
print("0.6 jacard = {}".format(d3 / (d3 +  d1 + d2)))
print("missing similarities fraction={}".format(d1/(d1+d3)))
print("extra similarities fraction={}".format(d2/(d2+d3)))