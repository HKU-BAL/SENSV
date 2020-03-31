import os
import logging
import argparse
import configparser
import pandas as pd
from threading import Thread
from pathlib import Path


def get_depth(samtools, chr, window, bam, output, logger):
    logger.info(f'Calculating mean depth of chr{chr}')
    cmd = (
        f'{samtools} depth -ar {chr} {bam} | '
        f'pypy3 fastmeandepthcnt.py - {window} > {output}'
    )
    os.system(cmd)

    logger.info(f'Done calculating mean depth for chr{chr}')


def depth_to_csv(output_name, output_dir, chr_list):
    df_files = (
        pd.read_csv(
            Path(output_dir) / f'depth_df_{chromo}.csv',
            names=["index", "start", "depth", "chr"],
            header=None
        ) for chromo in chr_list
    )
    df = pd.concat(df_files, ignore_index=True)
    df.to_csv(f'{output_dir}_{output_name}.csv', index=False, header=False)


def main(args=None):
    bam, output_path = args.bam, Path(args.output_path)
    # read config
    current_directory = Path(os.path.dirname(__file__)).resolve().parent.parent
    config_file_path = current_directory / "config.ini"
    config = configparser.ConfigParser()
    config.read(config_file_path)
    samtools = config['common']['samtools']

    # init logging stuff
    logger = logging.getLogger(name='log')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    name = os.path.basename(bam).split(".")[0]
    output_dir = f'chr_depth_{name}'
    os.makedirs(output_path / output_dir, exist_ok=True)

    chr_list = [str(x) for x in range(1, 23)] + ["X", "Y"]
    threads = [
        Thread(
            target=get_depth,
            kwargs={
                'samtools': samtools,
                'chr': chr,
                'window': '10000',
                'bam': bam,
                'output': f'{output_path / output_dir / f"depth_df_{chr}.csv"}',
                'logger': logger,
            }
        ) for chr in chr_list
    ]

    for t in threads:
        t.start()
    for t in threads:
        t.join()

    depth_to_csv(
        "all",
        output_path / output_dir,
        ["X", "Y"] + [str(x) for x in range(1, 23)]
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate mean depth from bam.')
    parser.add_argument('bam', help='input bam file.')
    parser.add_argument('output_path', help='output path.')
    args = parser.parse_args()

    main(args)
