import logging
import os
import argparse
import pandas as pd
from threading import Thread
import configparser

config = configparser.ConfigParser()
config.read("config.ini")

logger = logging.getLogger('log')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)


def getDepth(chr, window, bam, output):
    logger.info("Calculating mean depth of chr%s" % (chr))
    os.system(
        "samtools depth -ar %s %s | python fastmeandepthcnt.py - %s > %s" % (
            chr, bam, window, output))
    logger.info("Done calculating mean depth for chr%s" % (chr))


def depth_to_csv(output_name, output_dir, chr_list):
    df_files = (
    pd.read_csv("%s/depth_df_%s.csv" % (output_dir, chromo), names=["index", "start", "depth", "chr"], header=None) for
    chromo in chr_list)
    df = pd.concat(df_files, ignore_index=True)
    df.to_csv('%s_%s.csv' % (output_dir, output_name), index=False, header=False)


def main(args=None):
    parser = argparse.ArgumentParser(description='Calculate mean depth from bam.')
    parser.add_argument('bam', help='input bam file.')
    parser.add_argument('output_path', help='output path.')
    parser.add_argument('window',help='window size')
    args = parser.parse_args()
    name = os.path.basename(args.bam).split(".")[0]
    output_dir = "chr_depth_%s" % name
    os.makedirs(os.path.join(args.output_path, output_dir), exist_ok=True)
    threads = []
    chr_list = [str(x) for x in range(1, 23)] + ["X", "Y"]
    for chr in chr_list:
        t = Thread(target=getDepth, args=[chr, args.window, args.bam,
                                          "%s/depth_df_%s.csv" % (os.path.join(args.output_path, output_dir), chr)])
        threads.append(t)
        t.start()
    for t in threads:
        t.join()
    depth_to_csv("all", os.path.join(args.output_path, output_dir), ["X", "Y"] + [str(x) for x in range(1, 23)])
    # depth_to_csv("sex", output_dir, ["X", "Y"])


if __name__ == "__main__":
    main()
