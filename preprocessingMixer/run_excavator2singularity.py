#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from spython.main import Client
import os
import pandas as pd
import sys
import argparse
import shutil
import glob
import yaml


def guess_assembly(centromere_exca):
    ref = os.path.splitext(os.path.basename(centromere_exca))[0].split('_')[1]
    return ref    

def create_target_yaml(cen,fasta,bw,chromos,gap,tar,tmp):
    assem = guess_assembly(cen)
    fasta_dir = os.path.dirname(fasta)
    fasta_name = os.path.basename(fasta)
    bw_dir = os.path.dirname(bw)
    bw_name = os.path.basename(bw)
    chr_dir = os.path.dirname(chromos)
    chr_name = os.path.basename(chromos)
    cen_dir = os.path.dirname(cen)
    cen_name = os.path.basename(cen)
    gap_dir = os.path.dirname(gap)
    gap_name = os.path.basename(gap)    
    target_dir = os.path.dirname(tar)
    target_name = os.path.basename(tar)
    dir_name = os.path.splitext(os.path.basename(tar))[0]
    target_dict = {'Reference' : {'Assembly': assem, 'FASTA': '/efasta/'+fasta_name, 'BigWig':'/ebw/'+bw_name, 'Chromosomes': '/echr/'+chr_name,'Centromeres':'/ecen/'+cen_name,'Gaps': '/egap/'+gap_name},'Target' : {'Name':dir_name, 'BED':'/etarget/'+target_name, 'Window': int(50000)}}
    target_bind = fasta_dir+'/:/efasta/,'+bw_dir+'/:/ebw/,'+chr_dir+'/:/echr/,'+cen_dir+'/:/ecen/,'+gap_dir+'/:/egap/,'+target_dir+'/:/etarget/,'+tmp+'/:/output/'
    target_path = '/output/Target/'+target_dict['Reference']['Assembly']+'/'+target_dict['Target']['Name']+'/w_50000'  
    return target_dict,target_bind,target_path
    
def create_prepare_yaml(confile, bam_dir):    
    bam_bind = ''
    prepare_dict ={}
    startcount = confile.ID
    for i in range(len(startcount)): ##in case all bam are in different paths, bind each sample singularly
        bam_name = confile.bamPath[i]
        if i < len(startcount)-1:
           bam_bind += bam_dir+'/:/'+confile.ID[i]+'_dir/,'
        else:
           bam_bind += bam_dir+'/:/'+confile.ID[i]+'_dir/' 
        prepare_dict[confile.ID[i]] = '/'+confile.ID[i]+'_dir/'+bam_name
    return bam_bind,prepare_dict
    
def create_analysis_yaml(c,t):
    c_count = c.ID
    t_count = t.ID
    analysis_dict = {}
    for k in range(len(c_count)):
        analysis_dict["C"+str(k+1)] = c.ID[k]
    for j in range(len(t_count)):
        analysis_dict["T"+str(j+1)] = t.ID[j] 
    return analysis_dict 

def run_EXCA2(tyaml,pyaml,ayaml,b,Tpath,sd,ecpu, no_controls):
    tyaml_name = os.path.basename(tyaml)
    pyaml_name = os.path.basename(pyaml)
    ayaml_name = os.path.basename(ayaml)
    #Client.load(sd+"/"+"excavator2.sif")
    sif_path= os.path.join(sd, "excavator2.sif")
    Client.load(sif_path)
    print("Running EXCAVATOR2 TargetPerla")
    #Client.execute(['TargetPerla.pl','-v','-s','/output/'+tyaml_name,'-o','/output/Target'], bind=b, quiet=False)
    Client.run(sif_path,['TargetPerla.pl','-v','-s','/output/'+tyaml_name,'-o','/output/Target'],bind=b, quiet=False)
    print("Running EXCAVATOR2 DataPrepare")
    Client.run(sif_path, ['EXCAVATORDataPrepare.pl','-v','-s','/output/'+pyaml_name,'-o','/output/DataPrepare_w50k','-@',ecpu,'-t',Tpath], bind=b, quiet=False)
    #NOTE that mixer can use a control.RData, which can be provided
    if not no_controls:
        print("Running EXCAVATOR2 DataAnalysis")
        Client.run(sif_path, ['EXCAVATORDataAnalysis.pl','-v','-s','/output/'+ayaml_name,'-o','/output/DataAnalysis_w50k','-@',ecpu,'-t',Tpath, '-i','/output/DataPrepare_w50k', '-e','pooling'], bind=b, quiet=False)
    else:
        print("No control samples provided, skipping EXCAVATOR2 DataAnalysis")
    print("EXCAVATOR2 Analysis completed")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run EXCAVATOR2')
    parser.add_argument('-o', '--output_dir', metavar="", help="path to output dir")
    parser.add_argument('-b', '--bam_dir', help="path to bams directory.")
    parser.add_argument('-e', '--exp_name', default = "", metavar="", help="experiment or target name")
    parser.add_argument('-t', '--target', metavar="", help="target bed file - the same used for alignment and EXCAVATOR2")
    parser.add_argument('-cf', '--config', metavar="", help="miXer config file with XXXY samples")
    parser.add_argument('-r','--ref', metavar="", help="reference fasta file - the same used for alignment and EXCAVATOR2", default=False)
    #parser.add_argument('-r37','--ref37', metavar="", help="OPTIONAL: only with b37 build, add hg19 reference to run EXCAVATOR2", default=False)
    parser.add_argument('-m','--mapp', metavar="", help="EXCAVATOR2 bigwig mappability file")
    parser.add_argument('-g','--gap', metavar="", help="EXCAVATOR2 GAP regions file")
    parser.add_argument('-cm','--centromeres', metavar="", help="EXCAVATOR2 centromeres file")
    parser.add_argument('-ch','--chromosomes', metavar="", help="EXCAVATOR2 chromosomes file")
    parser.add_argument('-th', '--thread', metavar="", help="number of threads to use for EXCAVATOR2", type=str)
    parser.add_argument('-ed', '--excavator_dir', help="folder containing excavator2.sif file.")

    if len(sys.argv)==1:
        print()
        print("Usage: python", sys.argv[0]," --help")
        print()
        raise ValueError("No arguments provided")
    arguments = parser.parse_args()
    ########### DECLARE VARIABLES:
        
    ###tmp dir
    absP = os.path.abspath(arguments.output_dir)
    tmp_folder = os.path.join(absP,arguments.exp_name + "_excavator2_output/")
    if os.path.exists(tmp_folder) == False:
       os.makedirs(tmp_folder)
       print("EXCAVATOR2 temp folder created")
       
    ####inputs:
    samples_df = pd.read_table(arguments.config, sep="\t")
    samples_df.sampleType = samples_df.sampleType.str.replace(' ', '').str.lower()
    ##select from config files female controls
    samples_ctrl = samples_df[samples_df['sampleType'] == 'c'].reset_index(drop=True)
    no_controls = False
    if samples_ctrl.shape[0] == 0:
        no_controls = True

    samples_case = samples_df[samples_df['sampleType'] == 't'].reset_index(drop=True)
    if samples_case.empty:
        samples_case = samples_df[samples_df['sampleType'] == 'train'].reset_index(drop=True)[:1]
    ##script dir containing sif file
    #script_dir = os.path.realpath(os.path.dirname(__file__))


    #### create excavator2 yaml files:
    ##TargetPerla:
    #if arguments.ref37:
    #   target_dict,target_bind,target_path = create_target_yaml(arguments.centromeres,arguments.ref37,arguments.mapp,arguments.chromosomes,arguments.gap,arguments.target,tmp_folder)
    #else:
    target_dict,target_bind,target_path = create_target_yaml(arguments.centromeres,arguments.ref,arguments.mapp,arguments.chromosomes,arguments.gap,arguments.target,tmp_folder)
    target_yaml = os.path.join(tmp_folder,arguments.exp_name + '_target.yaml')
    with open(target_yaml, 'w') as file:
         documents = yaml.dump(target_dict, file, sort_keys=False)
    ##DataPrepare:
    bam_bind,prepare_dict = create_prepare_yaml(samples_df, arguments.bam_dir)
    prepare_yaml = os.path.join(tmp_folder,arguments.exp_name + '_dataprepare.yaml')
    with open(prepare_yaml, 'w') as file:
         documents = yaml.dump(prepare_dict, file, sort_keys=False)    
    ##DataAnalysis:
    analysis_dict = create_analysis_yaml(samples_ctrl,samples_case)
    analysis_yaml = os.path.join(tmp_folder,arguments.exp_name + '_analysis.yaml')
    with open(analysis_yaml, 'w') as file:
         documents = yaml.dump(analysis_dict, file, sort_keys=False)    
    ##concatenate all binds
    bind_all = target_bind + ','+ bam_bind
    
    ###run EXCAVATOR2
    run_EXCA2(target_yaml,prepare_yaml,analysis_yaml,bind_all,target_path,arguments.excavator_dir,arguments.thread, no_controls)

    ####rm tmp dir
    #shutil.rmtree(tmp_folder, ignore_errors=True) 


    

