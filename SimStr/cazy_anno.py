'''
Copyright {2020} Junyu Chen

# This file is part of StrainPanDA. 
#
# StrainPanDA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# StrainPanDA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with StrainPanDA.  If not, see <https://www.gnu.org/licenses/>.

'''

#updata log
#2022-05-20 by Yuxiang Tan
#add the checkpoint of the seqID format in seqFilter function to avoid error when processing WGS genes.
#add the parameter of faa, for WGS genes should use the faa file directly rather than running SeqkitTranslate
#update the way to get sample infor in the data processing section

import os
import re
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
#from shutil import copyfile
from itertools import repeat
from multiprocessing import Pool, freeze_support


def manifestGen(InDir, suffix):
    path = pd.DataFrame()
    for subdir, dirs, files in os.walk(InDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if file.endswith(suffix) and os.path.getsize(filePath) > 0:
                #print(filePath)
                sampleID = file.replace(suffix, "").replace("panphlan_", "")
                path = path.append({'SampleID':str(sampleID), "ffnList":str(filePath)}, ignore_index=True)
    return path

#filter seq len does not match to the annotation
def seqFilter(fileList, OutDir, suffix):
    check_warning=False
    for file in fileList:
        seqList = []
        seqErrList = []
        for seq in SeqIO.parse(file, "fasta"):
            loc = seq.id.split(":")[-1].split("-")
            if len(loc)==2:
                if abs(int(loc[1]) - int(loc[0])) + 1 == len(seq):
                    seqList.append(seq)
                else:
                    seqErrList.append(seq)
            else:
                check_warning=True
                seqList.append(seq)
        if check_warning:
            print("Warning: the seq id in input files are not in pangenome format. Should consider use the faa parameter instead.")
        
        #print(seqList)
        #print(os.path.join(OutDir, os.path.split(file)[1]))
        SeqIO.write(seqList, os.path.join(OutDir, os.path.split(file)[1]), "fasta")
        SeqIO.write(seqErrList, os.path.join(OutDir, os.path.split(file)[1].replace(suffix,"".join([suffix,"_err"]))), "fasta")

#run seqkit to translate seq in 6 frames
def RunSeqkitTranslate(inFile, OutDir, suffix):
    sample = os.path.split(inFile)[1].replace(suffix, ".faa")
    cmd = "seqkit translate --frame 6 --append-frame " + inFile + " -o " + os.path.join(OutDir, sample)
    subprocess.call(cmd, shell=True)

def RunSeqkitTranslateParallel(inFileList, OutDir, jobs, suffix):
    pool = Pool(processes = jobs)
    pool.starmap(RunSeqkitTranslate, zip(inFileList, repeat(OutDir), repeat(suffix)))
    pool.close()
    pool.join()
    pool.terminate()

## Run dbCAN
def RunDbCAN(faaFile, prefix, OutDir, databseDir, threads):
    cmd = "run_dbcan.py " + faaFile + " protein --out_dir " + os.path.join(OutDir, prefix) + " --db_dir " + databseDir + \
        " --dia_cpu " + str(threads) + " --hmm_cpu " + str(threads) + " --hotpep_cpu " + str(threads) + " --tf_cpu " + str(threads)
    subprocess.call(cmd, shell=True)
def RunDbCANParallel(faaList, prefixList, OutDir, databseDir, threads, jobs):
    pool = Pool(processes=jobs)
    pool.starmap(RunDbCAN, zip(faaList, prefixList, repeat(OutDir), repeat(databseDir), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

## parse dbCAN result
def parseDbcan_old(dbcanDir):
    filePathList = []
    for subdir, dirs, files in os.walk(dbcanDir):
        for file in files:
            filePath = os.path.join(subdir, file)
            if file.endswith(".txt") and os.path.getsize(os.path.join(subdir, file)) > 0:
                filePathList.append(filePath)
    frames = [ pd.read_table(f) for f in filePathList ]
    df_concat = pd.concat(frames)
    return df_concat

def parseDbcan(dbcanDir):
    filePathList = []
    frames = []
    dirs = os.listdir(dbcanDir)
    for dir in dirs:
        filePath = os.path.join(dbcanDir, dir, "overview.txt")
        if os.path.getsize(filePath) > 0:
            df_f=pd.read_table(filePath)
            df_f["SampleID"]=dir
            frames.append(df_f)
    df_concat = pd.concat(frames)
    return df_concat


# Function to clean the names 
def Clean_names(CAZy_name): 
    # Search for opening bracket in the name followed by 
    # any characters repeated any number of times 
    if re.search('\(.*', CAZy_name): 
        # Extract the position of beginning of pattern 
        pos = re.search('\(.*', CAZy_name).start() 
        # return the cleaned name 
        return CAZy_name[:pos] 
    else: 
        # if clean up needed return the same name 
        return CAZy_name 
# Function to clean hmm names 
def Clean_hmm_names(CAZy_name): 
    CAZy_name = str(CAZy_name).split("_")[0]
    return CAZy_name
def removeTail(string):
    pattern = "_frame=" + r'[0-9]'
    # Match all digits in the string and replace them by empty string
    mod_string = re.sub(pattern, '', string)
    return mod_string

parser = argparse.ArgumentParser(description='CAZy annotations')
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True,
                    help="the path of the reads")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True,
                    help="the output path of reads")
parser.add_argument('-d', '--database', dest='database', type=str,  required=False, default="/home/chenjunyu/Lab/Anno/database/dbCAN2_db",
                    help="the reference_reads path")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='8',
                    help="the number of threads run for a job")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='2',
                    help="the number of jobs run in parallel")
parser.add_argument('-s', '--suffix', dest='suffix', type=str, required=False, default="_centroids.ffn",
                    help="the suffix used to recognize the input files, if they are protein seqs, should add the -p option")
parser.add_argument('-p', '--faa', dest='faa', type=str, required=False, default="False",
                    help="the any input other than False will be consider True and skip the RunSeqkitTranslate step")
args = parser.parse_args()


InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
dbPath = os.path.abspath(args.database)
threads = int(args.threads)
jobs = int(args.jobs)
suffix = str(args.suffix)
faa = str(args.faa)

seqFilterDir = os.path.join(OutDir, "seq_clean")
if os.path.exists(seqFilterDir) == 0:
    os.makedirs(seqFilterDir, 0o777, True)
dbcanDir = os.path.join(OutDir, "dbcan")
if os.path.exists(dbcanDir) == 0:
    os.makedirs(dbcanDir, 0o777, True)
cazyDir = os.path.join(OutDir, "cazy")
if os.path.exists(cazyDir) == 0:
    os.makedirs(cazyDir, 0o777, True)

#######filter seq############
#suffix = "_centroids.ffn"
path = manifestGen(InDir, suffix)
path.to_csv(os.path.join(OutDir, "ffnPathTable_ori.csv"))
fnnList = path["ffnList"].tolist()
if faa=="False":
    seqFilter(fnnList, seqFilterDir, suffix)

###########translate seq#############
path1 = manifestGen(seqFilterDir, suffix)
path1.to_csv(os.path.join(OutDir, "ffnPathTable_filter.csv"))
ffnList = path1["ffnList"].tolist()
if faa=="False":
    RunSeqkitTranslateParallel(ffnList, seqFilterDir, jobs, suffix)

#############protein seq annotation #####################
if faa=="False":
    path2 = manifestGen(seqFilterDir, ".faa")
else:
    path2 = manifestGen(InDir, suffix)
faaList = path2["ffnList"].tolist()
prefixList = path2["SampleID"].tolist()
RunDbCANParallel(faaList, prefixList, dbcanDir, dbPath, threads, jobs)

#####data processing####################
df_concat = parseDbcan(dbcanDir)
df = df_concat.loc[df_concat["#ofTools"] >= 2]
df.loc[:, 'HMMER'] = df['HMMER'].apply(Clean_names)
df.loc[:, 'Hotpep'] = df['Hotpep'].apply(Clean_names)
df.loc[:, 'HMMER'] = df['HMMER'].apply(Clean_hmm_names)
df = df.reset_index()

#rename
for i in range(len(df)):
    if df.loc[i, "HMMER"] == "-":
        df.loc[i, "cazy"] = df.loc[i, "Hotpep"]
        if df.loc[i, "Hotpep"] == "-":
            df.loc[i, "cazy"] = df.loc[i, "DIAMOND"]
    else: df.loc[i, "cazy"] = df.loc[i, "HMMER"]

df["GeneID"] = df["Gene ID"].apply(removeTail)
df = df.drop_duplicates(subset=['GeneID'], keep='first')
#final out for all the samples together
final = pd.DataFrame()
final["SampleID"] = df["SampleID"]
final["GeneID"] = df["Gene ID"]
final["CAZy"] = df["cazy"]
final.to_csv(os.path.join(OutDir, "cazy_final_out.csv"), index=None)

#sep output for each sample
sampleList = list(df["SampleID"].unique())
for sample in sampleList:
    df1 = pd.DataFrame()
    df2 = df.loc[df["SampleID"] == sample]
    df1["GeneID"] = df2["GeneID"]
    df1["CAZy"] = df2["cazy"]
    df1.to_csv(os.path.join(cazyDir, sample + ".csv"), index=None)
