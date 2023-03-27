import sys
import os
import argparse
import subprocess
from pathlib import Path
#import configparser

#config = configparser.ConfigParser()
#config.read("config.ini")


def main(args=None):
	parser = argparse.ArgumentParser(description='get depth and normalize bam.')
	parser.add_argument('bam', nargs='?', help='input bam file.')
    #parser.add_argument('reference', nargs='?', help='input ref file.')
    #parser.add_argument('gender', nargs='?', help='m: male, f:female.')
	parser.add_argument('window', nargs='?', help='window size')
	parser.add_argument('output_path', nargs='?', help='output path.')
	parser.add_argument('random',nargs='?',help='randomly assign the bam to female and male')  
	args = parser.parse_args()
	os.chdir(sys.path[0])
	name = os.path.basename(args.bam).split(".",1)[0]
	male=[]
	female=[]
    #name = args.bam_num
	os.system("python3 depth.v1.py %s %s %s" % (args.bam, args.output_path,args.window))
	all_chr = subprocess.Popen("python3 normalize_sample_mask_somatic_norm_separately.py %s %s" % (os.path.join(args.output_path, "chr_depth_"+name+"_all.csv"), "mask_none.csv"), shell=True)
	all_chr.wait()
	#Y=os.path.join(args.output_path,"chr_depth_"+name,"depth_df_Y.csv")
	#count=0
	#f=open(Y,'r')
	#data=[]
	#for line in f:
	#	row=line.split(',')
	#	data.append(row)
	#	if(row[2]=="0.0"):
	#		count=count+1
	#if(float(count/len(data))>0.5):
	#	isMale=True
	#if(isMale):
	if(int(args.random)%2==0):
		male.append(os.path.join(args.output_path, "chr_depth_"+name+"_all_mask_somatic_norm.csv"))
	else:
		female.append(os.path.join(args.output_path, "chr_depth_"+name+"_all_mask_somatic_norm.csv"))
	path=Path(args.output_path)
	f=open(os.path.join(path.parent,"female.txt"),'a+')
	for i in range(0,len(female)):
		f.write(female[i]+"\n")
	m=open(os.path.join(path.parent,"male.txt"),'a+')
	for i in range(0,len(male)):
		m.write(male[i]+"\n")
    #all_chr = subprocess.Popen("python3 revised_with_input_thread.py %s %s %s %s" % (
    #os.path.join(args.output_path, "chr_depth_" + name + "_all_mask_somatic_norm.csv"), args.reference, args.gender, args.output_path),
    #                        shell=True)
    #all_chr.communicate()

if __name__ == "__main__":
	main()

