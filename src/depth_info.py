from argparse import ArgumentParser


class DepthInfoOptions:
    def __init__(self, in_csv, out_file):
        self.in_csv = in_csv
        self.out_file = out_file


class DepthInfo:
    """
    class DepthInfo is used to manipulate the results generated by depth distribution scripts
    """

    def __init__(self, options):
        self.in_csv = options.in_csv
        self.out_file = options.out_file

    # extract called regions from the resulting files generated by depth distribution scripts
    def extract_csv(self):
        start_line = "Match region:"
        endLine = ""

        depth_region = []
        with open(self.in_csv) as f:
            is_start = 0
            for line in f:
                if line.strip() == start_line:
                    is_start = 1
                    continue

                if is_start:
                    if line.strip() == endLine:
                        is_start = 0
                        break

                    arr = line.strip().split()
                    start, end, chrom = arr[0], arr[1], arr[2]
                    if chrom == '23':
                        chrom = 'X'
                    elif chrom == '24':
                        chrom = 'Y'

                    depth_region.append({'chrom': chrom, 'start': start, 'end': end})

        #self.depth_region = depth_region
        def toIntStr(text): return text if text.isdigit() else '23' if text == 'X' else '24'
        self.depth_region = sorted(depth_region, key=lambda k: toIntStr(k['chrom']).rjust(2) + str(k['start']).rjust(9))

    # generate output files
    def gen_output(self):
        with open(self.out_file, 'w') as f:
            for region in self.depth_region:
                print(region['chrom'], region['start'], region['end'], file=f)

    def run(self):
        self.extract_csv()
        self.gen_output()


if __name__ == "__main__":
    parser = ArgumentParser(description='run')
    parser.add_argument('-incsv', '--in_csv', help='input csv file', required=True)
    parser.add_argument('-outfile', '--out_file', help='output file', required=True)
    options = parser.parse_args()

    depth_info = DepthInfo(options)
    depth_info.run()
