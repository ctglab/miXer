#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#function that load mosdepth singularity image and compute the mean coverage over the autosomal target regions
#returns  ll mean coverages as single output in the given output directory

from spython.main import Client
import os
import pandas as pd
import sys
import argparse
import shutil
import glob

def fix_target(a):
    target_slice=a.iloc[:,:3]
    if any(c.isalpha() for c in str(target_slice.iloc[0,2])) == True:
       target_slice = target_slice.drop([0]).reset_index(drop=True)
    ###remove chrX,chrY,alt contigs
    target_bed=target_slice.loc[(target_slice.loc[:,0] != "chrX")]
    target_bed=target_bed.loc[(target_bed.loc[:,0] != "chrY")]
    target_bed=target_bed.loc[(target_bed.loc[:,0] != "chrM")]
    #Ensuring that all values in first column are string
    target_bed.iloc[:, 0] = target_bed.iloc[:, 0].astype(str)
    target_bed=target_bed[~target_bed.loc[:,0].str.contains('_')]
    if target_bed.empty:
       raise ValueError("Check that target file is properly formatted")
    else:
        return target_bed

def autosomal_cvg(tmpdir,sample_id,autosomaltar,bam,sd):
    input_dir = os.path.dirname(bam)
    bam_name = os.path.basename(bam)
    target_name = os.path.basename(autosomaltar)
    Client.load(os.path.join(sd, "mosdepth:v0.3.3.sif"))
    #Client.load(sd+"/"+"mosdepth:v0.3.3.sif")
    Client.execute(['mosdepth','-t','4','-x','-n','-b','/Tmp/'+target_name,'/Tmp/'+sample_id, '/inputdir/'+bam_name], bind=input_dir+'/:/inputdir/,'+tmpdir+'/:/Tmp/')
    os.remove(glob.glob(os.path.join(tmpdir,'*region.dist.txt'))[0])
    os.remove(glob.glob(os.path.join(tmpdir,'*global.dist.txt'))[0])
    os.remove(glob.glob(os.path.join(tmpdir,'*csi'))[0])
    os.remove(glob.glob(os.path.join(tmpdir,'*regions.bed.gz'))[0]) 
    
# def create_summary(outfile,cvg):
#     for cvgfile in cvg:
#         filename=os.path.basename(cvgfile).split(".")[0]
#         p=pd.read_csv(cvgfile, index_col=False, sep = "\t")
#         cvg = p.loc[p['chrom'] == "total_region", 'mean'].item()
#         outfile = outfile.append({'ID': filename, 'MeanCvg': cvg}, ignore_index=True)
#     return outfile

def create_summary(outfile, cvg):
    dfs = []
    for cvgfile in cvg:
        filename = os.path.basename(cvgfile).split(".")[0]
        p = pd.read_csv(cvgfile, index_col=False, sep="\t")
        cvg_value = p.loc[p['chrom'] == "total_region", 'mean'].item()
        df = pd.DataFrame({'ID': [filename], 'MeanCvg': [cvg_value]})
        dfs.append(df)
    return pd.concat([outfile] + dfs, ignore_index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute autosomal coverage for train and case samples using Mosdepth')
    parser.add_argument('-o', '--output_dir', metavar="", help="path to output dir")
    parser.add_argument('-e', '--exp_name', default = "", metavar="", help="experiment or target name")
    parser.add_argument('-t', '--target', metavar="", help="target bed file - the same used for alignment and EXCAVATOR2")
    parser.add_argument('-cf', '--config', metavar="", help="tab separated file with all samples autosomes coverage: ID meancvg", type=argparse.FileType('r'))
    parser.add_argument('-bd', '--bam_dir', help="directory containing BAM files reported in config file.")
    parser.add_argument('-md', '--mosdepth_dir', help="Directory of mosdepth.sif file.")
    
    if len(sys.argv)==1:
        print()
        print("Usage: python", sys.argv[0]," --help")
        print()
        raise ValueError("No arguments provided")
    arguments = parser.parse_args()

    ########### DECLARE VARIABLES:
    cvg_folder = None
    ###tmp dir
    tmp_folder = os.path.join(arguments.output_dir,"tmp_folder/")
    if os.path.exists(tmp_folder) == False:
       os.makedirs(tmp_folder)
       print("temp folder created")
    
    ####inputs:
    samples_df = pd.read_table(arguments.config, names= ['ID','bamPath','Gender','analysis'], sep="\t")
    samples_df["bamPath"] = samples_df['bamPath'].apply(lambda x: os.path.join(arguments.bam_dir, x))
    samples_df.analysis = samples_df.analysis.str.replace(' ', '').str.lower()
    samples_df.ID = samples_df.ID.astype(str)
    ##select from config files sample to process for both training and calling
    samples = samples_df[samples_df['analysis'].str.contains(r'train')]
    target_df=pd.read_csv(arguments.target, sep="\t", index_col=False, header=None)
    target_autosomal = fix_target(target_df)
    target_out = os.path.join(tmp_folder,'target_tmp.bed')
    target_autosomal.to_csv(target_out, index=None, header=None, sep="\t")
    #script_dir = os.path.realpath(os.path.dirname(__file__))
    ## start mosdepth
    if not samples.size == 0:
       start = samples['bamPath']
       for i in range(len(start)):
           #autosomal_cvg(tmp_folder,samples['ID'][i],target_out,samples['bamPath'][i],script_dir)        
           autosomal_cvg(tmp_folder,samples['ID'][i],target_out,samples['bamPath'][i],arguments.mosdepth_dir)
           
       ##collect cvgs
       cvg_folder = os.path.join(arguments.output_dir,arguments.exp_name + "_autosomal_cvg.txt")
       mosdepth_list = list(glob.glob(os.path.join(tmp_folder,'*summary.txt')))
       outDf = pd.DataFrame(columns = ['ID','MeanCvg'])
       create_summary(outDf,mosdepth_list).to_csv(cvg_folder,index=None, header=None, sep="\t")
    
    ###rm tmp dir
    shutil.rmtree(tmp_folder, ignore_errors=True) 
    if cvg_folder is not None:
        if os.path.getsize(cvg_folder) > 0:
            print("Coverage calculation completed")
        else:
           raise ValueError("Coverage file is empty, check that mosdepth summary files were created")



  
    

