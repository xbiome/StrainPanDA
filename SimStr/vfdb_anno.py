'''
Copyright {2020} Yuxiang Tan

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

import os
import re
import argparse
import subprocess
import pandas as pd
from Bio import SeqIO
from itertools import repeat
from multiprocessing import Pool, freeze_support


####logic and steps: 
#1. use manifestGen and seqFilter to find all target files by suffix, into the PathTable.csv;
#2. For DNA seq, use RunSeqkitTranslateParallel to generate protein seqs, or for proteins skipping this step;
#3. RunDiamondParallel to run diamond onto target database
#4. for pangenome, use parse_6orf_blastp, select the highest pident for each gene; for protein seq, use parseBlastp.
#The output is in the final_out folder with one file per sample.
#There are two columns in the file, the first column is the gene/gf ID, and the second column is the ID in the database.
#An All_samples_features_out.csv file is generated as the merged of all samples with all features.

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

#run diamond
def RunDiamond(fasta, prefix, dbDir, blastpOut, threads):
    out_format = '6 qseqid sseqid pident qcovhsp length mismatch gapopen qlen qstart qend sstart send evalue bitscore'
    cmd = "diamond blastp --db " + dbDir + " --query " + fasta + " --out " + os.path.join(blastpOut, prefix + "_blasp.tsv") + " --evalue 1e-05 --outfmt " + out_format + " --max-target-seqs 1" + " --threads " + str(threads)
    subprocess.call(cmd, shell=True)

def RunDiamondParallel(fileList, prefixList, dbDir, blastpOut, threads, jobs):
    pool = Pool(processes = jobs)
    pool.starmap(RunDiamond, zip(fileList, prefixList, repeat(dbDir), repeat(blastpOut), repeat(threads)))
    pool.close()
    pool.join()
    pool.terminate()

#for Seqkit results
def removeTail(string):
    pattern = "_frame=" + r'[0-9]'
    # Match all digits in the string and replace them by empty string
    mod_string = re.sub(pattern, '', string)
    return mod_string


#parse functions
#This is for 6 orf from pangenome
def parse_6orf_blastp(blastpOut, pident, qcov,blastp_finalout):
    #out = pd.DataFrame(columns = ['qseqid', 'sseqid', 'pident', 'qcovhsp', 'length', 'mismatch','gapopen', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue','bitscore', 're_qseqid', 'SampleID', 'GeneFamily'])
    out = pd.DataFrame()
    for file in os.listdir(blastpOut):
        if file.endswith("_blasp.tsv") and os.path.getsize(os.path.join(blastpOut, file)) > 0:
            filePath = os.path.join(blastpOut, file)
            SampleID =file.replace("_blasp.tsv", "")
            df = pd.read_table(filePath, header=None)
            df.columns = ["qseqid", "sseqid", "pident", "qcovhsp", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            df50 = df.loc[(df["pident"] > pident) & (df["qcovhsp"] > qcov)]
            #6 orf processing
            df50.loc[:, "re_qseqid"] = df["qseqid"].apply(removeTail)
            df50 = df50.sort_values(by = ["pident"], ascending=False)
            df50_dedup = df50.drop_duplicates(subset=["re_qseqid"], keep='first')
            if len(df50_dedup) > 0:
                df50_dedup.loc[:, "SampleID"] = SampleID
                filePath_out=os.path.join(blastp_finalout, file.replace("_blasp.tsv", "_anno.csv"))
                out = out.append(df50_dedup)
                final = pd.DataFrame()
                final["GeneID"] = df50_dedup["re_qseqid"]
                final["FuncID"] = df50_dedup["sseqid"]
                final.to_csv(filePath_out, index=None)
                
    return out

#This is for faa input, and hard coded identity
def parseBlastp(blastpOut, pident, qcov,blastp_finalout):
    df1 = pd.DataFrame()
    for file in os.listdir(blastpOut):
        if file.endswith("_blasp.tsv") and os.path.getsize(os.path.join(blastpOut, file)) > 0:
            filePath = os.path.join(blastpOut, file)
            SampleID =file.replace("_blasp.tsv", "")
            df = pd.read_table(filePath, header=None)
            df.columns = ["qseqid", "sseqid", "pident", "qcovhsp", "length", "mismatch", "gapopen", "qlen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
            df50 = df.loc[(df["pident"] > pident) & (df["qcovhsp"] > qcov)]
            if len(df50) > 0:
                df50.loc[:, "SampleID"] = SampleID
                filePath_out=os.path.join(blastp_finalout, file.replace("_blasp.tsv", "_anno.csv"))
                df1 = df1.append(df50, ignore_index=True)
                final = pd.DataFrame()
                final["GeneID"] = df50["qseqid"]
                final["FuncID"] = df50["sseqid"]
                final.to_csv(filePath_out, index=None)
    return df1


parser = argparse.ArgumentParser(description='Welcome Run blast')
parser.add_argument('-i', '--input', dest='InDir', type=str, required=True, help="the path of the genome file")
parser.add_argument('-o', '--output', dest='OutDir', type=str, required=True, help="the output path of result")
parser.add_argument('-j', '--jobs', dest='jobs', type=str,  required=False, default='4', help="the number of jobs run in parallel")
parser.add_argument('-db', '--database', dest='dbDir', type=str,  required=True, help="the path of the database")
parser.add_argument('-pi', '--pident', dest='pident', type=str,  required=False, default='30', help="the cutoff value of identity")
parser.add_argument('-cov', '--cov', dest='cov', type=str,  required=False, default='30', help="the cutoff value of coverage")
parser.add_argument('-x', '--extension', dest='suffix', type=str,  required=False, default='_centroids.ffn', help="the extension of input file")
parser.add_argument('-p', '--faa', dest='faa', type=str, required=False, default="False",
                    help="the any input other than False will be consider True and skip the RunSeqkitTranslate step")
parser.add_argument('-t', '--threads', dest='threads', type=str, required=False, default='8',
                    help="the number of threads run for a job")
args = parser.parse_args()

InDir = os.path.abspath(args.InDir)
OutDir = os.path.abspath(args.OutDir)
dbDir = os.path.abspath(args.dbDir)
pident = int(args.pident)
qcov = int(args.cov)
suffix = str(args.suffix)
jobs = int(args.jobs)
threads = int(args.threads)
faa = str(args.faa)

if os.path.exists(OutDir) == 0:
    os.makedirs(OutDir, 0o777, True)

seqFilterDir = os.path.join(OutDir, "seq_clean")

if os.path.exists(seqFilterDir) == 0:
    os.makedirs(seqFilterDir, 0o777, True)

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

#run diamond
if faa=="False":
    path2 = manifestGen(seqFilterDir, ".faa")
else:
    path2 = manifestGen(InDir, suffix)

path2.to_csv(OutDir + "/PathTable.csv", index=False, encoding = "utf-8")
faaList = path2["ffnList"].tolist()
prefixList = path2["SampleID"].tolist()

blastpOut = os.path.join(OutDir, "blastp")
if os.path.exists(blastpOut) == 0:
    os.makedirs(blastpOut, 0o777, True)

RunDiamondParallel(faaList, prefixList, dbDir, blastpOut, threads, jobs)

blastp_finalout = os.path.join(OutDir, "final_out")
if os.path.exists(blastp_finalout) == 0:
    os.makedirs(blastp_finalout, 0o777, True)   


#now the target is to generate relation between gene ID or gf ID with function, more advance merging is not consider.
if faa=="False":
    df = parse_6orf_blastp(blastpOut, pident, qcov,blastp_finalout)
else:
    df = parseBlastp(blastpOut, pident, qcov,blastp_finalout)
df.to_csv(os.path.join(OutDir, "All_samples_features_out.csv"), index=None)




