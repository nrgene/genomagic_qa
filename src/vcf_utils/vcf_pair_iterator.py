import sys


def parse_vcf_line(line):
    parts = line.split('\t')
    chr1 = int(parts[0])
    start1 = int(parts[1])
    assert parts[7][:7] == 'HS;END='
    end1 = int(parts[7][7:])
    return [chr1, start1, end1]


def get_first_non_header_line(my_file):
    first_line='#'
    while first_line[0] == '#':
        first_line = my_file.readline().rstrip()
    return first_line


def advance_one_file_and_get_line(progeny_file):
    line1 = progeny_file.readline().rstrip()
    assert line1
    return line1


def iterate_two_vcf(first_file, second_file, line1, line2, vcf_lines_intersection_method):

    while line1 or line2:
        assert line1
        [chr1, start1, end1] = parse_vcf_line(line1)
        assert line2
        [chr2, start2, end2] = parse_vcf_line(line2)
        if chr1 < chr2 or (chr1 == chr2 and start2 >= end1):
            line1 = advance_one_file_and_get_line(first_file)
        elif chr1 > chr2 or (chr1 == chr2 and start1 >= end2):
            line2 = advance_one_file_and_get_line(second_file)
        else:
            #at this stage, we know that we have intersection between the two vcf lines
            assert chr1 == chr2
            assert start1 < end2
            assert start2 < end1
            vcf_lines_intersection_method(line1, line2)
            if end1 == end2:
                line1 = first_file.readline().rstrip()
                line2 = second_file.readline().rstrip()
            elif end1 < end2:
                line1 = advance_one_file_and_get_line(first_file)
            else:
                assert end2 < end1
                line2 = advance_one_file_and_get_line(second_file)

