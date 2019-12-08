
def register_data_version_command(host, dv_name):
    url = 'http://{}/genomagic-api/v1/data-versions/{}'.format(host, dv_name)
    return 'curl -X PUT \'{}\''.format(url)


def add_samples_group_command(samples_list, group_name, host, dv_name):
    samples_string = ",".join(["\"{}\"".format(s) for s in samples_list])
    sample_string_json = "{}\"sampleIds\" : [{}]{}".format("{", samples_string, "}")
    url = "http://{}/genomagic-api/v1/sample-groups/{}?dataVersion={}".format(host, group_name, dv_name)
    my_command = " curl -X PUT \'{}\' - H \'Content-Type: application/json\' -d \'{}\'".format(url, sample_string_json)
    return my_command