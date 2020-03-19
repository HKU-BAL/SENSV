#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import argparse

COLUMN_TYPES = [str, int, int] # chrom, position, depth

parser = argparse.ArgumentParser(
    description='Calculate mean depth from samtools depth output.')
parser.add_argument('depth', type=argparse.FileType('r'), nargs='?',
                    help='samtools depth output.')
parser.add_argument('window', type=int, help='window size')
args = parser.parse_args()

region_start = None
depth_sum = 0
depth_count = 0
output_count = 0

for line in args.depth:
    chrom, pos, depth = map(lambda conv, v: conv(v), COLUMN_TYPES,
                            line.strip().split('\t'))
    if region_start is None:
        region_start = pos
    elif pos >= region_start + args.window and depth_count > 0 and args.window > 0:
        print(output_count, str(pos), str(depth_sum/depth_count), chrom, sep=',')
        output_count += 1
        region_start = pos
        depth_sum = 0
        depth_count = 0
    depth_sum += depth
    depth_count += 1

if depth_count > 0:
    print(output_count, str(region_start + depth_count), str(depth_sum/depth_count), chrom, sep=',')
