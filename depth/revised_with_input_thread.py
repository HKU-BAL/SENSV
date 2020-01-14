#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:30:06 2019

@author: lavenderyu
"""
import argparse
from sys import maxsize
import numpy as np
import pandas as pd
import pickle
import csv
import math
import scipy.stats
import os
from pathlib import Path
import sys
from multiprocessing import Pool
parser = argparse.ArgumentParser(description='Calculate likelihood')

parser.add_argument('sample_norm_file', type=str, help='sample norm file path for likelihood calculation')
parser.add_argument('reference', help='input ref file.')
parser.add_argument('gender',type=str,help='patient gender')
parser.add_argument('output_dir', help='output directory.')
args = parser.parse_args()
ref = []
regions = []
sample = []
sample_file = args.sample_norm_file
reference_len = 0
sample_len = 0
with open(args.reference, "r") as f:
    has_header = csv.Sniffer().has_header(f.readline())
    f.seek(0)
    readRef = csv.reader(f, delimiter=',')

    if has_header:
        next(readRef)
    j = 0
    for line in readRef:
        ref.append([i if i == "X" or i == "Y" or i == "NA" or i == "" else float(i) for i in line])
        j += 1
    reference_len = j

    assert reference_len != 0, "the reference is empty"

with open(sample_file, "rU") as f:
    has_header = csv.Sniffer().has_header(f.readline())
    f.seek(0)
    readSample = csv.reader(f, delimiter=',')

    if has_header:
        next(readSample)
    j = 0
    for line in readSample:
        sample.append([i if i == "X" or i == "Y" or i == "NA" or i == "" else float(i) for i in line])
        j += 1
    sample_len = j
    print("sample length: %i" % sample_len)
    print("reference length: %i" % reference_len)
    assert sample_len == reference_len, "the sample should have the same length as reference"


def local_max_subarray(df, size, min_len=500000):
    candidates = []
    max_ending_here = 0 #df.at[0, 'diflognorm']
    temp_start_index = temp_end_index = None
    prev_n_candidates = 0

    # find maximum subarray by greedy
    curPos = 1

    start2row = {}
    for row in range(size):
        start2row[df.at[row, 'start']] = row

    while curPos < size - 1:
        curDif = df.at[curPos, 'diflognorm']
        curPos += 1
        
        if curDif > (max_ending_here + curDif):
            # new start of summation

            if len(candidates) > prev_n_candidates:
                # rewind current position to end of last peak if we have a new peak
                curPos = start2row[int(candidates[-1][-1])] 
                temp_start_index = temp_end_index = df.at[curPos, 'start']
                max_ending_here = 0
                prev_n_candidates = len(candidates)
                continue
            else:
                temp_start_index = temp_end_index = df.at[curPos, 'start']
                max_ending_here = curDif
            continue

        # extend the existing summation
        temp_end_index = df.at[curPos if curPos != size else size - 1, 'start']
        max_ending_here = max_ending_here + curDif
        if (temp_start_index == None): temp_start_index = df.at[0, 'start']
        

        # skip updating candidates array if subarray length is not long enough
        if temp_end_index == None or temp_start_index == None or temp_end_index - temp_start_index < min_len: continue

        # if candidates is empty or the last item is a brand new subarray => append candidates
        if len(candidates) == 0 or candidates[-1][1] != temp_start_index:
            candidates.append([max_ending_here, temp_start_index, temp_end_index])
        # else if the last candidate is not local maximum
        elif candidates[-1][1] == temp_start_index and max_ending_here >= candidates[-1][0]:
            candidates[-1] = [max_ending_here, temp_start_index, temp_end_index]

    #print (candidates)
    return candidates


# create a dataFrame for dif of log likelihood between two models
def getDifByChr(chr_num):
    diflikhood_half = pd.DataFrame(columns=['start', 'diflognorm'])
    diflikhood_dup=pd.DataFrame(columns=['start','diflognorm'])
    for row in range(len(ref)):
        if (ref[row][22] == chr_num and chr_num != 23 and chr_num!=24):

            xs=[np.float32(ref[row][1]),np.float32(ref[row][2]),np.float32(ref[row][3]),np.float32(ref[row][4]),np.float32(ref[row][5]),np.float32(ref[row][6]),np.float32(ref[row][7]),np.float32(ref[row][8]),np.float32(ref[row][9]),np.float32(ref[row][10]),np.float32(ref[row][11]),np.float32(ref[row][12]),np.float32(ref[row][13]),np.float32(ref[row][14]),np.float32(ref[row][15]),np.float32(ref[row][16]),np.float32(ref[row][17]),np.float32(ref[row][18]),np.float32(ref[row][19]),np.float32(ref[row][20])]
            xs_w=scipy.stats.mstats.winsorize(xs,limits=[0,0.05])
            mean=sum(xs_w)/len(xs_w)
            if(mean==0): continue;
            std=np.std(xs_w)
            pos = ref[row][21]
            x = np.float32(sample[row][1])
            norm = scipy.stats.norm(mean, std).pdf(x)
            halfnorm = scipy.stats.norm(mean /2, std).pdf(x)
            dupnorm=scipy.stats.norm(mean * 1.5,std).pdf(x)
            if (norm == 1):
                meanlog=0
            else:
                if (norm==0):
                   norm=sys.float_info.min * sys.float_info.epsilon
                meanlog=math.log(norm)
            if (halfnorm == 1):
                halfmeanlog=0
            else:
                if(halfnorm==0):
                    halfnorm=sys.float_info.min * sys.float_info.epsilon
                halfmeanlog = math.log(halfnorm)
            if (dupnorm == 1):
                dupmeanlog=0
            else:
                if (dupnorm == 0):
                    dupnorm = sys.float_info.min * sys.float_info.epsilon
                dupmeanlog=math.log(dupnorm)

            dif_half = halfmeanlog - meanlog
            dif_dup=dupmeanlog-meanlog
            df_half= pd.DataFrame([[pos, dif_half]], columns=['start', 'diflognorm'])
            df_dup = pd.DataFrame([[pos, dif_dup]],columns=['start','diflognorm'])
            diflikhood_half = diflikhood_half.append(df_half, ignore_index=True)
            diflikhood_dup=diflikhood_dup.append(df_dup,ignore_index=True)

        elif(ref[row][22] ==chr_num and (chr_num==23 or chr_num==24)):#26
            #if (depth_3 =='y' and args.gender =='m'):
            #    xs=[np.float32(ref[row][1]),np.float32(ref[row][2]),np.float32(ref[row][3]),np.float32(ref[row][9]),np.float32(ref[row][11]),np.float32(ref[row][14]),np.float32(ref[row][15]),np.float32(ref[row][25]),np.float32(ref[row][26]),np.float32(ref[row][35]),np.float32(ref[row][39]),np.float32(ref[row][42]),np.float32(ref[row][46]),np.float32(ref[row][49])]
            #elif(depth_3=='n' and args.gender=='m'):
            #    xs=[np.float32(ref[row][4]),np.float32(ref[row][7]),np.float32(ref[row][17]),np.float32(ref[row][20]),np.float32(ref[row][23])]
            #elif(depth_3=='y' and args.gender=='f'):
            #    xs=[np.float32(ref[row][48]),np.float32(ref[row][47]),np.float32(ref[row][45]),np.float32(ref[row][10]),np.float32(ref[row][44]),np.float32(ref[row][43]),np.float32(ref[row][16]),np.float32(ref[row][41]),np.float32(ref[row][40]),np.float32(ref[row][38]),np.float32(ref[row][37]),np.float32(ref[row][36])]                
            #else:
            #    xs=[np.float32(ref[row][5]),np.float32(ref[row][6]),np.float32(ref[row][8]),np.float32(ref[row][12]),np.float32(ref[row][13]),np.float32(ref[row][18]),np.float32(ref[row][19]),np.float32(ref[row][21]),np.float32(ref[row][22]),np.float32(ref[row][24])]               
            if(args.gender=='f'):
                xs=[np.float32(ref[row][1]),np.float32(ref[row][2]),np.float32(ref[row][3]),np.float32(ref[row][4]),np.float32(ref[row][7]),np.float32(ref[row][9]),np.float32(ref[row][8]),np.float32(ref[row][6]),np.float32(ref[row][5]),np.float32(ref[row][10])]
            elif(args.gender=='m'):
                xs=[np.float32(ref[row][17]),np.float32(ref[row][20]),np.float32(ref[row][15]),np.float32(ref[row][14]),np.float32(ref[row][16]),np.float32(ref[row][13]),np.float32(ref[row][18]),np.float32(ref[row][19]),np.float32(ref[row][12]),np.float32(ref[row][11])]  
            xs_w=scipy.stats.mstats.winsorize(xs,limits=[0,0.05])                
            mean=sum(xs_w)/len(xs_w)     
            if(mean==0):
                    continue;
            std=np.std(xs_w)
            pos=ref[row][21]
            x=sample[row][1]           
            norm=scipy.stats.norm(mean, std).pdf(x)
            halfnorm=scipy.stats.norm(mean/2, std).pdf(x)
            dupnorm=scipy.stats.norm(mean * 1.5,std).pdf(x)
            if (norm == 1):
                meanlog=0
            else:
                if (norm==0):
                   norm=sys.float_info.min * sys.float_info.epsilon
                meanlog=math.log(norm)
            if (halfnorm == 1):
                halfmeanlog=0
            else:
                if(halfnorm==0):
                    halfnorm=sys.float_info.min * sys.float_info.epsilon
                halfmeanlog = math.log(halfnorm)
            if (dupnorm == 1):
                dupmeanlog=0
            else:
                if (dupnorm == 0):
                    dupnorm = sys.float_info.min * sys.float_info.epsilon
                dupmeanlog=math.log(dupnorm)

            dif_half = halfmeanlog - meanlog
            dif_dup=dupmeanlog-meanlog
            df_half= pd.DataFrame([[pos, dif_half]], columns=['start', 'diflognorm'])
            df_dup = pd.DataFrame([[pos, dif_dup]],columns=['start','diflognorm'])
            diflikhood_half = diflikhood_half.append(df_half, ignore_index=True)
            diflikhood_dup=diflikhood_dup.append(df_dup,ignore_index=True)
    return diflikhood_half,diflikhood_dup
    
def process(chrm):
    results_del = pd.DataFrame(columns=['chr', 'start', 'end','sum of diflihood', 'del'])
    results_dup = pd.DataFrame(columns=['chr', 'start', 'end','sum of diflihood', 'dup'])

    difByChr_half,difByChr_dup = getDifByChr(chrm)
    regions_half = local_max_subarray(difByChr_half, len(difByChr_half.index), 500000)

    regions_dup=local_max_subarray(difByChr_dup,len(difByChr_dup.index),500000)

    for i in range(0, len(regions_half)):
        likelihood=regions_half[i][0]
        likelihood_chi_square=scipy.stats.chi2.sf(2*likelihood,1)
        if(likelihood_chi_square <0.05):
          df = pd.DataFrame([[chrm, regions_half[i][1], regions_half[i][2],likelihood, "DEL"]], columns=['chr', 'start', 'end', 'sum of diflihood','del'])
          results_del = results_del.append(df, ignore_index=True)

    for i in range(0, len(regions_dup)):
        likelihood=regions_dup[i][0]
        likelihood_chi_square=scipy.stats.chi2.sf(2*likelihood,1)
        if(likelihood_chi_square <0.05):
          df = pd.DataFrame([[chrm, regions_dup[i][1], regions_dup[i][2],likelihood, "DUP"]], columns=['chr', 'start', 'end', 'sum of diflihood','dup'])
          results_dup = results_dup.append(df, ignore_index=True)

    return {'results_del':results_del, 'results_dup':results_dup}

def main():
    # os.makedirs(args.output_dir, exist_ok=True)
    # drive code
    results_del = pd.DataFrame(columns=['chr', 'start', 'end','sum of diflihood', 'del'])
    results_dup=pd.DataFrame(columns=['chr','start','end','sum of diflihood','dup'])
    end=-1
    if args.gender=='f':
        end=24
    elif args.gender=='m':
        end=25

    chrm_list = range(1, end)
    pool = Pool(len(chrm_list))
    results = pool.map(process, chrm_list)
    pool.close()
    pool.join()


    for result in results:
         results_del = results_del.append(result['results_del'], ignore_index=True)
         results_dup = results_dup.append(result['results_dup'], ignore_index=True)

    output_prefix = ""
    if args.output_dir=="./" or args.output_dir=="/" or args.output_dir==".":
        output_prefix = os.path.basename(args.sample_norm_file).split(".")[0]
    else:
        if args.output_dir[-1] == "/":
            output_prefix = args.output_dir[:-1].rsplit("/", 1)[1]
        elif "/" in args.output_dir:
            output_prefix = args.output_dir.rsplit("/", 1)[1]
        else:
            output_prefix=args.output_dir
    results_del.to_csv(os.path.join(os.path.dirname(args.output_dir), output_prefix + ".depth"), index=None, header=False, sep='\t')
    with open(os.path.join(os.path.dirname(args.output_dir), output_prefix + ".depth"),'a') as f:
        results_dup.to_csv(f,index=None,header=False, sep='\t')

if __name__ == "__main__":
    main()

