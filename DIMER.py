#!/usr/bin/env python
# Original Programmer: jianan.lin@jax.org
# Version: 1.0.1
# Publishing date: October 6, 2017
# Required packages: pybedtools, numpy, scipy, math, and pandas.
import sys
import multiprocessing
import pybedtools
from pybedtools.featurefuncs import *
import os
import os.path
from subprocess import *
import time
import numpy as np
import math
import scipy
from scipy.stats import norm
import pandas as pd
import argparse
from argparse import RawTextHelpFormatter


def show_value(s):
    if sys.version_info.major == 2:
        if isinstance(s, unicode):
            return str(s)
    return s

def featuretype_filter(feature, featuretype):
    return feature[2] == featuretype

def subset_featuretypes(featuretype,bedtoolobj):
    result = bedtoolobj.filter(featuretype_filter, featuretype).saveas()
    return pybedtools.BedTool(result.fn)


def gene_level_bed(bedname,outname):
    ifile=open(bedname)
    ofile=open(outname,'w')
    for line in ifile:
        [chrom,start,end,name,score,strand]=line.strip("\r\n").split("\t")
        new_chr=chrom+"___"+name
        ofile.write("\t".join([new_chr,start,end,name+strand,score,strand])+"\n")
    ofile.close()
    ifile.close()

def modify_ID(ifilename,ofilename):
    ifile=open(ifilename)
    ofile=open(ofilename,'w')
    for line in ifile:
        [new_chr,start,end,name]=line.strip("\r\n").split("\t")
        chrom=new_chr.split("___")[0]
        ofile.write("\t".join([chrom,start,end,name[:-1],".",name[-1]])+"\n")
    ofile.close()
    ifile.close()

def get_read_counts(bamfilename):
    bam=pybedtools.BedTool(bamfilename)
    countfile=bamfilename[:-4]+".count.txt"
    bed.coverage(bam.fn,counts=True,s=True).saveas(countfile)

def combine_countfile(bamfile_list,ofilename):
    ls=[filename[:-4]+".count.txt" for filename in bamfile_list]
    d={}
    for filename in ls:
        ifile=open(filename)
        for line in ifile:
            [chrom,start,end,name,score,strand,count]=line.strip("\r\n").split("\t")
            try:
                d['___'.join([chrom,start,end,name,strand])].append(count)
            except KeyError:
                d['___'.join([chrom,start,end,name,strand])] = [count]
        ifile.close()
    ofile=open(ofilename,'w')
    ofile.write("\t".join(["chr","start","end","name","strand"]+[x[:-4] for x in ls])+"\n")
    for key in d:
        ofile.write('\t'.join(key.split("___")+d[key])+'\n')
    ofile.close()
    return ls


def get_splited_trans(gtf_name,window_size,fofilename,bamfile_list):
    ## Step 1: get exons ##
    ofilename=gtf_name[:-4]+".win"+str(window_size)+".transcriptome.bed"
    if os.path.isfile(ofilename)==False:
        print "Build Transcriptome from gtf file ..."
        g=pybedtools.BedTool(gtf_name)
        exons=subset_featuretypes('exon',g)
        exonbedfilename=gtf_name[:-4]+".bed"
        exons.each(gff2bed).saveas(exonbedfilename)
        gene_exon_bed_name=gtf_name[:-4]+".gene.bed"
        gene_level_bed(exonbedfilename,gene_exon_bed_name)
        ge=pybedtools.BedTool(gene_exon_bed_name)
        ## Step 2: Split exons ##
        print "Split exons ..."
        split_tranname=gtf_name[:-4]+".split.gene.bed"
        gemerge=ge.sort().merge(c=4,o="distinct")
        gemerge.window_maker(b=gemerge.fn,w=window_size,i="src").saveas(split_tranname)
        #ofilename=gtf_name[:-4]+".transcriptome.bed"
        modify_ID(split_tranname,ofilename)
        cmd="rm {0}\nrm {1}\nrm {2}".format(exonbedfilename,gene_exon_bed_name,split_tranname)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
    ## Step 3: Get counts from bam files ... ##
    print "Calculating coverage ...."
    if os.path.isfile(fofilename)==False:
        global bed
        bed=pybedtools.BedTool(ofilename)
        #pool=multiprocessing.Pool()
        #pool.map(get_read_counts,bamfile_list)
        map(get_read_counts,bamfile_list)
        print "Combine read count files ..."
        ls=combine_countfile(bamfile_list,fofilename)
        for filename in ls:
            cmd="rm {0}".format(filename)
            Popen(cmd,shell=True,stdout=PIPE).communicate()

def get_meth_lv(IP_ls,Input_ls):
    if len(IP_ls)==len(Input_ls):
        meth_ls=[float(x)/float(x+y+(1e-12)) for x,y in zip(IP_ls,Input_ls)]
        return meth_ls
    else:
        return False

def cqtest(X,Y):
    Result = [0,0]
    (n1,p1)=X.shape
    (n2,p2)=Y.shape
    if not p1==p2:
        print "Dimensions do not match"
        return "Dimension Error"
    p=p1
    if (n1<=2) or (n2<=2):
        print "Minimum sample size required for both groups is 3."
        return "Dimension Error"
    if np.sum(np.isnan(X)) + np.sum(np.isnan(Y)) > 0:
        Result=[0,1]
    elif (np.sum(X**2)==0) and (np.sum(Y**2)==0):
        Result=[0,1]
    elif (np.sum(X**2)!=0) and (np.sum(Y**2)==0):
        t=0
        d=0
        q=range(n1)
        for i in q:
            for j in [q[x] for x in q if x!=i]:
                t+=np.sum(X[i,]*X[j,])
                if n1>3:
                    Xl=np.delete(X,[i,j],axis=0)
                    Amean=np.mean(Xl,axis=0)
                else:
                    Amean=np.delete(X,[i,j],axis=0)
                d+=np.sum(X[j,]*(X[i,]-Amean))*np.sum(X[i,]*(X[j,]-Amean))
        if d!=0:
            a=(t/float(n1*(n1-1)))/math.sqrt(2*d/float((n1*(n1-1))**2))
            b=1-norm.cdf(a,loc=0,scale=1)
            Result=[a,b]
        elif d==0:
            Result=[0,1]
    elif (np.sum(X**2)==0) and (np.sum(Y**2)!=0):
        t=0
        d=0
        q=range(n2)
        for i in q:
            for j in [q[x] for x in q if x!=i]:
                t+=np.sum(Y[i,]*Y[j,])
                if n2>3:
                    Yl=np.delete(Y,[i,j],axis=0)
                    Amean=np.mean(Yl,axis=0)
                else:
                    Amean=np.delete(Y,[i,j],axis=0)
                d+=np.sum(Y[j,]*(Y[i,]-Amean))*np.sum(Y[i,]*(Y[j,]-Amean))
        if d!=0:
            a=(t/float(n2*(n2-1)))/math.sqrt(2*d/float((n2*(n2-1))**2))
            b=1-norm.cdf(a,loc=0,scale=1)
            Result=[a,b]
        elif d==0:
            Result=[0,1]
    else:
        t1=0
        t2=0
        t3=0
        d1=0
        d2=0
        d3=0
        q1=range(n1)
        q2=range(n2)
        for i in q1:
            for j in [q1[x] for x in q1 if x!=i]:
                t1+=np.sum(X[i,]*X[j,])
                if n1>3:
                    Xl=np.delete(X,[i,j],axis=0)
                    Amean=np.mean(Xl,axis=0)
                else:
                    Amean=np.delete(X,[i,j],axis=0)
                d1+=np.sum(X[j,]*(X[i,]-Amean))*np.sum(X[i,]*(X[j,]-Amean))
        for i in q2:
            for j in [q2[x] for x in q2 if x!=i]:
                t2+=np.sum(Y[i,]*Y[j,])
                if n2>3:
                    Yl=np.delete(Y,[i,j],axis=0)
                    Amean=np.mean(Yl,axis=0)
                else:
                    Amean=np.delete(Y,[i,j],axis=0)
                d2+=np.sum(Y[j,]*(Y[i,]-Amean))*np.sum(Y[i,]*(Y[j,]-Amean))
        for i in q1:
            for j in q2:
                t3+=np.sum(X[i,]*Y[j,])
                Xl=np.delete(X,(i),axis=0)
                Yl=np.delete(Y,(j),axis=0)
                Amean1=np.mean(Xl,axis=0)
                Amean2=np.mean(Yl,axis=0)
                d3+=np.sum(Y[j,]*(X[i,]-Amean1)) * np.sum(X[i,]*(Y[j,]-Amean2))
        Numer=(t1/(n1*(n1-1))+t2/(n2*(n2-1))-2*t3/(n1*n2))
        Denomin=np.sqrt(2*d1/((n1*(n1-1))**2) + 2*d2/((n2*(n2-1))**2) + 4*d3/((n1*n2)**2))
        if Denomin!=0:
            a=Numer/float(Denomin)
            b=1-norm.cdf(a,loc=0,scale=1)
            Result=[a,b]
        elif Denomin==0:
            Result=[0,1]
    return Result


def meth_samples(fofilename,index_filename,indir):
    ifile=open(index_filename)
    d_index={}
    group_set={}
    next(ifile)
    for line in ifile:
        tmp=line.strip("\r\n").split()
        if len(tmp)>4:
            print "Please make sure there is no space exisiting in the path to your bam files.\nBad directory:\n"+line
            sys.exit()
        try:
            group_set[tmp[0]].append(tmp[0]+"_"+tmp[2])
        except KeyError:
            group_set[tmp[0]]=[tmp[0]+"_"+tmp[2]]
        if tmp[1]=="IP":
            try:
                d_index[tmp[0]+"_"+tmp[2]][0]=indir+tmp[3][:-4]+".count"
            except KeyError:
                d_index[tmp[0]+"_"+tmp[2]]=[indir+tmp[3][:-4]+".count",'']
        elif tmp[1]=="Input":
            try:
                d_index[tmp[0]+"_"+tmp[2]][1]=indir+tmp[3][:-4]+".count"
            except KeyError:
                d_index[tmp[0]+"_"+tmp[2]]=['',indir+tmp[3][:-4]+".count"]
    ifile.close()
    group_list=list(group_set)
    d={}
    d_meth={}
    d_meth["chr"]=[]
    d_meth["start"]=[]
    d_meth["end"]=[]
    ifile=open(fofilename)
    fl=ifile.readline().strip("\r\n").split("\t")
    dgene_ind={}
    for col in fl[5:]:
        d[col]=[]
    for i,line in enumerate(ifile):
        tmp=line.strip("\r\n").split("\t")
        d_meth["chr"].append(tmp[0])
        d_meth["start"].append(tmp[1])
        d_meth["end"].append(tmp[2])
        ggname=tmp[3]
        try:
            dgene_ind[ggname].append(i)
        except KeyError:
            dgene_ind[ggname]=[i]
        for col,j in zip(fl[5:],range(5,len(tmp))):
            d[col].append(int(tmp[j]))
    ifile.close()
    print "Finished reading count file"
    '''
    if os.path.isfile(fofilename[:-4]+"_meth.txt")==False:
        print "Generating meth table..."
        ofile=open(fofilename[:-4]+"_meth.txt",'w')
        ofile.write("chr\tstart\tend\tname\tstrand\t")
        for key in d_index:
            [IPname,Inputname]=d_index[key]
            d_meth[key]=get_meth_lv(d[IPname],d[Inputname])
            ofile.write(key+'\t')
        ofile.write("\n")
        for i,item in enumerate(d['___'.join(fl[:5])]):
            itemls=item.split('___')
            gene_list.append(itemls[3])
            ofile.write('\t'.join(map(str,itemls+[d_meth[key][i] for key in d_meth]))+'\n')
        ofile.close()
    '''
    for key in d_index:
        [IPname,Inputname]=d_index[key]
        d_meth[key]=get_meth_lv(d[IPname],d[Inputname])
        omethfile=fofilename[:-10]+"_meth_lv_"+key.split("/")[-1]+".bedgraph"
        df=pd.DataFrame(d_meth)
        df_cols=["chr","start","end",key]
        df.to_csv(omethfile,header=False,columns=df_cols,sep="\t",index=False)
        cmd="cat {0} | sort -k1,1 -k2,2n > {1}".format(omethfile,omethfile+"1")
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        cmd="cat {0} > {1}".format(omethfile+"1",omethfile)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
        cmd="rm {0}".format(omethfile+"1")
        Popen(cmd,shell=True,stdout=PIPE).communicate()
    print "Finished calculating the meth level."
    #otestfile=open(fofilename[:-4]+"_gene_test2.txt",'w')
    d_test={}
    for gene_name in dgene_ind:
        #print gene_name
        col_index=dgene_ind[gene_name]
        gary_signal=[]
        gary_control=[]
        for key in d_meth:
            ml=[d_meth[key][x] for x in col_index]
            if key in group_set[group_list[0]]:
                gary_signal.append(ml)
            elif key in group_set[group_list[1]]:
                gary_control.append(ml)
        gary_signal_use=np.array(gary_signal)
        gary_control_use=np.array(gary_control)
        d_test[gene_name]=cqtest(gary_signal_use,gary_control_use)
    print "Finished calculating test statistics and p-values."
    otestfile=open(fofilename[:-4]+"_gene_cqtest.txt",'w')
    for key in d_test:
        otestfile.write('\t'.join(map(str,[key]+d_test[key]))+'\n')
    #return d_test


def run_meth(gtf_name,window_size,fofilename,index_filename,indir):
    bamfile_list=[]
    ifile=open(index_filename)
    next(ifile)
    for line in ifile:
        tmp=line.strip("\r\n").split()
        if len(tmp)>4:
            print "Please make sure there is no spaces exisiting in the path to your bam files.\nBad directory:\n"+line
            sys.exit()
        bamfile_list.append(indir+tmp[3])
    ifile.close()
    get_splited_trans(gtf_name,window_size,fofilename,bamfile_list)
    meth_samples(fofilename,index_filename,indir)
    ofilename=gtf_name[:-4]+".win"+str(window_size)+".transcriptome.bed"
    if os.path.isfile(ofilename)==True:
        cmd="rm {0}".format(ofilename)
        Popen(cmd,shell=True,stdout=PIPE).communicate()

def gargparser():
    argparser = argparse.ArgumentParser(description="DIMER: DIfferential MEthylation of RNA\nVersion: 1.0", formatter_class=RawTextHelpFormatter)
    argparser.add_argument("-g","--gtf",dest = "gtffile", type = str, required = True, help = "input gtf file")
    argparser.add_argument("-i","--info",dest = "indexfile", type = str,required = True, help = "input PATH file. Please refer to the readme file for the format of this PATH file.")
    argparser.add_argument("-o","--outdir",dest="outdirectory", type = str,required=False,default="./", help = "output directory to save all ouput files. Default is the current folder")
    argparser.add_argument("-w","--winsize",dest = "winsize", type = int ,required = False, default=50, help = "window size used to create bins. Default is 50.")
    argparser.add_argument("-p","--pathin",dest = "indir", type = str,required = False, help = "The directory of the input bam files. If you have the valid directory information already in the file of -i option, you don't need to use this option. If you only specify the filename in the file of -i option, you need to provide this option for DIMER to find the input files.")
    return(argparser)

def mainfun():
    argsp=gargparser()
    args = argsp.parse_args()
    gtf_name=args.gtffile
    window_size=args.winsize
    if not os.path.exists(args.outdirectory):
        cmd="mkdir {0}".format(args.outdirectory)
        Popen(cmd,shell=True,stdout=PIPE).communicate()
    fofilename=args.outdirectory+"/win_"+str(window_size)+"_count.txt"
    index_filename=args.indexfile
    if args.indir is not None:
        indir=args.indir
        if not indir.endswith("/"):
            indir+="/"
    else:
        indir=""
    start_time = time.time()
    run_meth(gtf_name,window_size,fofilename,index_filename,indir)
    print "Finished."
    print("--- The running time is %s seconds ---" % (time.time() - start_time))


mainfun()


