import sys
genomagic_qa_repo_path = '/home/nrgsupport2/genomagic_qa_ariel/clients/genomagic_qa'
sys.path.append('{}/src'.format(genomagic_qa_repo_path))
import similarities_comparison
import api_server_utils
dv1 = "mon_maize_v2_3"
dv2 = "mon_maize_v3_4"
selected_chr = "1"
#samples = ["91inh2__ver110", "94ink1a__ver110"]
samples = ["3iih6__ver110", "53dwq1__ver100", "91inh2__ver110", "94ink1a__ver110", "biqa347__ver100", "gs2145__ver100", "lh172__ver100", "neto127__ver100", "sr1502__ver100", "yuca151_p11_dh1__ver110"]
api_server = "api.genomagic.local"
hap_sim1 = api_server_utils.get_raw_similarities_between_multiple_sampels(api_server, dv1, samples, selected_chr)
hap_sim2 = api_server_utils.get_raw_similarities_between_multiple_sampels(api_server, dv2, samples, selected_chr)
hap_sim1
[a1,a2,a12] = similarities_comparison.compute_similarity_match_for_multiple_sampels_in_df(hap_sim1, hap_sim2, samples)
print("{}, {}, {}".format(a1,a2,a12))
a_or_b = a1 + a2 + a12
print((1.0*a12) / a_or_b)
