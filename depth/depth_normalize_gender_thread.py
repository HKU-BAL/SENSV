import sys
import os
import argparse
import subprocess
import configparser

config = configparser.ConfigParser()
config.read("../config.ini")


def main(args=None):
    parser = argparse.ArgumentParser(description='get depth and normalize bam.')
    parser.add_argument('bam', nargs='?', help='input bam file.')
    parser.add_argument('reference', nargs='?', help='input ref file.')
    parser.add_argument('gender', nargs='?', help='m: male, f:female.')
    parser.add_argument('output_path', nargs='?', help='output path.')

    args = parser.parse_args()
    os.chdir(sys.path[0])
    name = os.path.basename(args.bam).rsplit(".",1)[0]
    os.system("python3 depth.thread.py %s %s" % (args.bam, args.output_path))
    all_chr = subprocess.Popen("python3 normalize_sample_mask_somatic_norm_separately.py %s %s" % (
    os.path.join(args.output_path, "chr_depth_" + name + "_all.csv"), config["depth"]["mask_table"]), shell=True)
    all_chr.communicate()
    all_chr = subprocess.Popen("python3 revised_with_input_thread.py %s %s %s %s" % (
    os.path.join(args.output_path, "chr_depth_" + name + "_all_mask_somatic_norm.csv"), args.reference, args.gender, args.output_path),
                               shell=True)
    all_chr.communicate()

if __name__ == "__main__":
    main()
