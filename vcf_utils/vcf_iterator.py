class VcfIterator:

    @staticmethod
    def get_last_header_and_first_non_header(input_file):
        line = input_file.readline().rstrip()
        assert line[0] == '#'
        header_line = ''
        while line[0] == '#':
            header_line = line
            line = input_file.readline().rstrip()
        return [header_line.split('\t'), line.split('\t')]

    # instance attribute
    def __init__(self, filename):
        self.vcf_file = open(filename, 'r')
        self.has_next_line = True
        [self.columns_header, self.next_line] = VcfIterator.get_last_header_and_first_non_header(self.vcf_file)

    def get_next_line_splitted(self):
        if not self.has_next():
            return None
        #assert self.has_next()
        curr_line = self.next_line
        my_next_line = self.vcf_file.readline().rstrip()
        if not my_next_line:
            self.has_next_line = False
            self.vcf_file.close()
        else:
            self.next_line = my_next_line.split('\t')
        return curr_line

    def has_next(self):
        return self.has_next_line

    def get_headers(self):
        return self.columns_header
