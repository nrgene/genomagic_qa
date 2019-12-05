def create_sample_string(sample_name, file1, file2):
    return "\"{}\":[\"{}\",\"{}\"]".format(sample_name, file1, file2)


def create_multiple_samples_string(arr):
    my_str = ''
    for x in arr:
        assert len(x) == 3
        sample_string = create_sample_string(x[0], x[1], x[2])
        my_str = '{},\n{}'.format(my_str, sample_string)
    return my_str[2:]


def create_sample_to_files_string(progeny_files):
    multiple_samples_string = create_multiple_samples_string(progeny_files)
    return "\"sampleToFiles\":{}\n{}\n{}".format("{", multiple_samples_string, "}")


def create_population_of_parents(parent1, parent2, progeny_files):
    parentalLines = "\"parentalLines\":[\"{}\",\"{}\"]".format(parent1, parent2)
    sampleToFiles = create_sample_to_files_string(progeny_files)
    return "{}{},{}{}".format("{", parentalLines, sampleToFiles, "}")


def create_populations_string_without_pedigree_file(parents_list, progeny_files_list):
    n = len(parents_list)
    assert len(progeny_files_list) == n
    all_parents_list = []
    for i in range(n):
        parents_pair = parents_list[i]
        assert len(parents_pair) == 2
        parents_str = create_population_of_parents(parents_pair[0], parents_pair[1], progeny_files_list[i])
        all_parents_list.append(parents_str)
    all_parents_string = ",\n".join(all_parents_list)
    return "\"populations\":[\n{}]".format(all_parents_string)


def create_population_string(pedigreeFileId, parents_list, progeny_files_list):
    populations_string_without_pedigree_file = create_populations_string_without_pedigree_file(parents_list, progeny_files_list)
    return "\"pedigreeFileId\":\"{}\",\n{}".format(pedigreeFileId, populations_string_without_pedigree_file)


def get_full_json_str(prev_dv, new_dv, pop_str):
    dv_str = "\"sourceDataVersion\":\"{}\",\n\"targetDataVersion\":\"{}\",".format(prev_dv, new_dv)
    return "{}\n{}\n{}\n{}".format("{", dv_str, pop_str, "}")


#gets - data version names,
# parents list - list of pairs - each pair is parent 1 parent2
# progeny_files_list - list in the length of the prev list, each element is a list of triplets, each is sample name, fastq1, fastq2
def create_parental_json(prev_dv, new_dv, pedigree_file_id, parents_list, progeny_files_list):
    pop_str = create_population_string(pedigree_file_id, parents_list, progeny_files_list)
    return get_full_json_str(prev_dv, new_dv, pop_str)


#gets - data version names, list of triplets, each is sample name, fastq1, fastq2
def create_diversity_json(prev_dv, new_dv, progeny_files):
    pop_str = create_sample_to_files_string(progeny_files)
    return get_full_json_str(prev_dv, new_dv, pop_str)