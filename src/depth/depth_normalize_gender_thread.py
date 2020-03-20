import os
import sys
import shlex
import argparse
import subprocess
import configparser
from pathlib import Path


def main(args=None):
    bam, reference, gender, output_path, ref_ver, nprocs = \
        args.bam, args.reference, args.gender, args.output_path, args.ref_ver, args.nprocs

    # current directory = repo directory (SENSV folder)
    current_directory = Path(os.path.dirname(__file__)).resolve().parent.parent
    mask_table_absolute_path = str(current_directory / 'data' / 'depth' / f'GRCh{ref_ver}.mask.csv')

    name = Path(bam).stem

    cmd = f'python depth.thread.py {bam} {output_path}'
    # print(f'cmd: #{cmd}#')
    os.system(cmd)

    cmd = (
        f'python normalize_sample_mask_somatic_norm_separately.py '
        f'{Path(output_path) / f"chr_depth_{name}_all.csv"} '
        f'{mask_table_absolute_path}'
    )
    # print(f'cmd: #{cmd}#')
    process = subprocess.Popen(shlex.split(cmd), stdout=sys.stdout, stderr=sys.stderr)
    process.communicate()

    cmd = (
        f'python revised_with_input_thread.py '
        f'{Path(output_path) / f"chr_depth_{name}_all_mask_somatic_norm.csv"} '
        f'{reference} {gender} {output_path} {nprocs}'
    )
    # print(f'cmd: #{cmd}#')
    process = subprocess.Popen(shlex.split(cmd), stdout=sys.stdout, stderr=sys.stderr)
    process.communicate()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='get depth and normalize bam.')
    parser.add_argument('bam', nargs='?', help='input bam file.')
    parser.add_argument('reference', nargs='?', help='input ref file.')
    parser.add_argument('gender', nargs='?', help='m: male, f:female.')
    parser.add_argument('output_path', nargs='?', help='output path.')
    parser.add_argument('ref_ver', nargs='?', help='ref_ver')
    parser.add_argument('nprocs', nargs='?', help='max # of processes to run.')
    args = parser.parse_args()

    main(args)
