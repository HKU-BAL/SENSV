import os
import argparse
import subprocess
parser=argparse.ArgumentParser(description='generate all depth from all bam')
parser.add_argument('window',help='window size')
parser.add_argument('output_dir',help='output directory for all the depth files')
parser.add_argument('bam_files',help='the file path that containing all the absolute path of the bam files')
args=parser.parse_args()

bamfiles=[]

bams=open(args.bam_files,'r').readlines()
for bam in bams:
	bamfiles.append(bam.rstrip())
#for root,dirs,files in os.walk(args.root,topdown=True):
#	for file in files:
#		if file.endswith('.bam'):
#			bamfiles.append(os.path.join(root,file))

i=0
while i<len(bamfiles):
	output_dir=args.output_dir
	name = os.path.basename(bamfiles[i]).split(".")[0]
	output=output_dir+name
	print("python gen_depth_samtools.py"+" "+bamfiles[i]+" "+args.window+" "+output)
	os.system("python gen_depth_samtools.py %s %s %s %s" % (bamfiles[i],args.window,output,i))
	print("samtools "+bamfiles[i]+" done")
	i=i+1

depth_file=os.path.join(output,"chr_depth_"+name+"_all.csv")
#poschr_path="poschr.csv"
poschr_path=os.path.join(args.output_dir,"poschr.csv")
#cmd="echo 'pos,chr' > poschr_path"
#os.system(cmd)
f=open(poschr_path,'w+')
f.write('pos,chr'+'\n')
o=open(depth_file)
for line in o:
	line=line.strip()
	row=line.split(',')
	if(row[3]=="X"):
		add=row[1]+","+"23"+"\n"
	elif(row[3]=="Y"):
		add=row[1]+","+"24"+"\n"
	else:
		add=row[1]+","+row[3]+"\n"
	f.write(add)
f.close()
o.close()

#command="cat "+depth_file+"|awk -F',' '{if($4==\"X\") {print $2\",\"\"23\"} else if($4==\"Y\") {print $2\",\"\"24\"} else {print $2\",\"$4}}' >> "+poschr_path+"temp"
#os.system(command)
#os.system("cat "+poschr_path+" "+poschr_path+"temp > "+os.path.join(args.output_dir,"poschr.csv"))
#print("cat "+poschr_path+" "+poschr_path+"temp > "+os.path.join(args.output_dir,"poschr.csv"))
print("poschr.csv created")
#print("start generating reference, could take 10 minutes")
print("python generate_reference.py "+os.path.join(output_dir))
ref_gen=subprocess.Popen("python generate_reference.py %s" % (os.path.join(args.output_dir)),shell=True)
#os.system("python generate_reference.py "+os.path.join(args.output_dir))
ref_gen.wait()
print("all is done")
