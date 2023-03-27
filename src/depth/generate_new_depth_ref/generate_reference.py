import os
import pandas
import csv
import numpy
import scipy.stats

import argparse
print("generating reference")
parser = argparse.ArgumentParser(description='generate reference')
parser.add_argument('root',help='root for all the csv files')
#parser.add_argument('output',help='csv output file')
args=parser.parse_args()

pos_chr=[]

depthfiles_f=[]
depthfiles_m=[]
female_files=os.path.join(args.root,"female.txt")
male_files=os.path.join(args.root,"male.txt")
f=open(female_files,'r')
for line in f:
	line=line.strip()
	depthfiles_f.append(line)
m=open(male_files,'r')
for line in m:
	line=line.strip()
	depthfiles_m.append(line)
df1=pandas.read_csv(os.path.join(args.root,"poschr.csv"))
#df1=pandas.read_csv("/nas8/hjyu/CNV/test_chun/depth_test_2/poschr.csv")
#print(df1)
print(len(df1.index))
#print(os.path.join(args.root,"poschr.csv"))
#for root,dirs,files in os.walk(args.root,topdown=True):
#  for file in files:
#    if file.endswith('norm.csv'):
#      depthfiles.append(os.path.join(root,file))
def genMatrix(depthfiles):
	refAll =[]
	i=0
	while i<len(depthfiles):
		depthSC=[]
		with open(depthfiles[i], "r") as f:
			has_header = csv.Sniffer().has_header(f.readline())
			f.seek(0)
			readFile = csv.reader(f, delimiter=',')
			if has_header:
				next(readFile)
			for line in readFile:
				depthSC.append(float(line[1]))
		refAll.append(depthSC)
		i=i+1
	refAllT=numpy.transpose(refAll)
	return refAllT

def generate_ref(female,male):
	#print("hihi")
	ref=[]
	ref_f=genMatrix(female)
	ref_m=genMatrix(male)
	#print(len(ref_f))
#	print(ref_m)
	if(len(ref_f)!=len(ref_m)):
		print("female ref is not the same with male, did you generate them based on different window sizes?")
		return ref
	if(len(ref_f)!=len(ref_m) or len(df1.index)!=len(ref_f) or len(df1.index)!=len(ref_m)):
		return ref
	m=0
	xs=[]
	ref=pandas.DataFrame(columns=['mean_n','std_n','mean_f','std_f','mean_m','std_m'])
	#print("hihi")
	for row in range(len(ref_f)):
		#print(df1.iloc[row][1])  
		if(df1.iloc[row][1]!=23 and df1.iloc[row][1]!=24):
			xs=ref_f[row]+ref_m[row]
			xs_w=scipy.stats.mstats.winsorize(xs,limits=[0,0.05])
			mean=sum(xs_w)/len(xs_w)
			std=numpy.std(xs_w)
			df=pandas.DataFrame([[mean, std,"","","",""]], columns=['mean_n','std_n','mean_f','std_f','mean_m','std_m'])		
		elif(df1.iloc[row][1]==23 or df1.iloc[row][1]==24):
			#print("haha")
			xs_f=ref_f[row]
			xs_fw=scipy.stats.mstats.winsorize(xs_f,limits=[0,0.05])
			mean_f=sum(xs_fw)/len(xs_fw)
			std_f=numpy.std(xs_fw)
			xs_m=ref_m[row]
			xs_mw=scipy.stats.mstats.winsorize(xs_m,limits=[0,0.05])
			mean_m=sum(xs_mw)/len(xs_mw)
			std_m=numpy.std(xs_mw)
			df=pandas.DataFrame([["", "",mean_f,std_f,mean_m,std_m]], columns=['mean_n','std_n','mean_f','std_f','mean_m','std_m'])
		#elif(df1.iloc[row][1]==24):
		#	xs=ref_m[row]
		#	xs_w=scipy.stats.mstats.winsorize(xs,limits=[0,0.05])
		#	mean=sum(xs_w)/len(xs_w)
		#	std=numpy.std(xs_w)
		#	df=pandas.DataFrame([["","","","",mean,std]], columns=['mean_n','std_n','mean_f','std_f','mean_m','std_m'])
		ref=ref.append(df,ignore_index=True)
	#print("done")
	return ref
#print(len(ref.index))
ref=generate_ref(depthfiles_f,depthfiles_m)
#print("done")
#print(len(ref.index))
#numpy.savetxt(arg.root+"temp.csv",ref,delimiter=",")
#df1=pandas.read_csv("/nas8/hjyu/CNV/simulation_ref_time.csv")
#os.system("cat "+depthfiles[0]+"|awk -F',' '{if($4==\"X\") {print $2\",\"\"23\"} else if($4==\"Y\") {print $2\",\"\"24\"} else {print $2\",\"$4}}' >> poschr.csv")
#df1=pandas.read_csv(os.path.join(args.root,"poschr.txt"))
#if(len(ref_f.index)==len(ref_m.index)):
#	ref=pandas.concat([ref_f,ref_m],axis=1)
if(len(ref.index)==len(df1.index)):
  	pandas.concat([ref,df1],axis=1).to_csv(os.path.join(args.root,"reference.csv"),index=False,na_rep='N/A')	
	#ref.to_csv(os.path.join(args.root,"reference.csv"),index=False,na_rep='N/A')  
