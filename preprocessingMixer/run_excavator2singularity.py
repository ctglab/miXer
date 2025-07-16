#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import shutil
import glob
import json
import logging

# setup logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler('excavator2.log'),
        logging.StreamHandler(sys.stdout)
    ])

try:
    import pandas as pd
except ModuleNotFoundError:
    logging.info(f"Pandas not found. Installing in the environment..")
    subprocess.run("python3 -m pip install pandas", shell=True)
    import pandas as pd
try:
    import yaml
except ModuleNotFoundError:
    logging.info(f"PyYAML not found. Installing in the environment..")
    subprocess.run("python3 -m pip install PyYAML", shell=True)
    import yaml

def guess_assembly(centromere_exca):
    ref = os.path.splitext(os.path.basename(centromere_exca))[0].split('_')[1]
    return ref    

<<<<<<< HEAD
def create_target_yaml(config: dict, tmp: str, window_size=50000) -> str:
    """
    the structure of the yaml from excavator2 is the following
    # Reference genome
    Reference:
        Assembly: hg38 # assembly name
        FASTA: .test/ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta # path to reference genome FASTA file
        BigWig: .test/ref/k100.Bismap.MultiTrackMappability.bw # path to reference genome BigWig file
        Chromosomes: .test/ref/ChromosomeCoordinate_hg38.txt # path to reference genome chromosome coordinates file
        Centromeres: .test/ref/CentromerePosition_hg38.txt # path to reference genome centromere positions file
        Gaps: .test/ref/GapHg38.UCSC.txt # path to reference genome gap positions file

    # Target regions
    Target:
        Name: SureSelectV7 # target name
        BED: .test/data/SureSelect_Human_All_Exon_V7.hg38.bed # path to target BED file
        Window: 30000 # sliding window size (bp)
    """
    target_dict = {'Reference': {}, 'Target': {}}
    assembly = guess_assembly(config['centro'])
    target_name = os.path.basename(config['target'])
    for yaml_key, json_key in zip(
        ['Assembly', 'FASTA', 'BigWig', 'Chromosomes', 'Centromeres', 'Gaps'],
        ['_', 'ref', 'map', 'chrom', 'centro', 'gap']
    ):
        if yaml_key == 'Assembly':
            target_dict['Reference'][yaml_key] = assembly
=======
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
        bam_name = confile.bamName[i]
        if i < len(startcount)-1:
           bam_bind += bam_dir+'/:/'+confile.ID[i]+'_dir/,'
>>>>>>> main
        else:
            target_dict['Reference'][yaml_key] = config[json_key]
    target_dict['Target']['Name'] = target_name
    target_dict['Target']['BED'] = config['target']
    target_dict['Target']['Window'] = window_size   
    target_path = os.path.join(
        tmp, 'output', 'Target', target_dict['Reference']['Assembly'], target_dict['Target']['Name'], f'w_{window_size}')
    target_yaml = os.path.join(
        tmp, f"{config['exp_id']}_target.yaml"
    )
    logging.info(f"Creating target YAML for excavator2 under {target_yaml}")
    with open(target_yaml, 'w') as f:
        yaml.dump(target_dict, f, sort_keys=False)
    return target_path, target_yaml

def create_prepare_yaml(config: dict, tmp: str) -> str:
    """
    The structure of the file is the following
    Test1:  .test/data/PL11_T_shallow_0.001.bam
    Control1: .test/data/PL11_N_shallow_0.01.bam
    """
    samples_info = pd.read_csv(config["sample_list"],
        sep="\t", dtype=str, skiprows=1, header=None, usecols=[0,1], names=[
        "ID", "bamPath"])
    prepare_yaml = os.path.join(tmp, f"{config['exp_id']}_dataprepare.yaml")
    dest_dict = samples_info.set_index('ID')['bamPath'].to_dict()
    logging.info(f"Creating prepare YAML for excavator2 under {prepare_yaml}")
    with open(prepare_yaml, 'w') as f:
        yaml.dump(dest_dict, f, sort_keys=False)
    return prepare_yaml

def create_analysis_yaml(config: dict, tmp: str) -> str:
    """
    This is just 
    T1: Test1
    C1: Control1
    """
    samples_info = pd.read_csv(config["sample_list"],
        sep="\t", dtype=str, skiprows=1, usecols=[0,3], header=None, names=[
        "ID", "sampleType"])
    analysis_yaml = os.path.join(tmp, f"{config['exp_id']}_analysis.yaml")
    samples_info['encoded_name'] = samples_info\
        .groupby('sampleType').cumcount()\
        .add(1).astype(str).radd(samples_info['sampleType']) # this creates T1, T2, .., C1
    samples_info = samples_info.loc[:, ['ID', 'encoded_name']]
    logging.info(f"Creating analysis YAML for excavator2 under {analysis_yaml}")
    with open(analysis_yaml, 'w') as f:
        mapping = samples_info.set_index("encoded_name")["ID"].to_dict()
        yaml.dump(mapping, f, sort_keys=False)
    return analysis_yaml    

def run_EXCA2(tyaml: str, pyaml: str, ayaml: str, Tpath: str, ecpu: int, no_controls):
    # tyaml_name = os.path.basename(tyaml)
    # pyaml_name = os.path.basename(pyaml)
    # ayaml_name = os.path.basename(ayaml)
    logging.info("Running EXCAVATOR2 TargetPerla")
    subprocess.run(['TargetPerla.pl','-v','-s', tyaml ,'-o', os.path.join(
        os.path.dirname(os.path.abspath(tyaml)), 'output/Target')])
    logging.info("Running EXCAVATOR2 DataPrepare")
    logging.info('Running ' + ' '.join(['EXCAVATORDataPrepare.pl','-v','-s', pyaml, '-o', os.path.join(
        os.path.dirname(os.path.abspath(pyaml)), 'output/DataPrepare_w50k'),'-@', str(ecpu),'-t',Tpath]))
    subprocess.run(f'cat {pyaml}', shell=True)
    subprocess.run(['EXCAVATORDataPrepare.pl','-v','-s', pyaml, '-o', os.path.join(
        os.path.dirname(os.path.abspath(pyaml)), 'output/DataPrepare_w50k'),'-@', str(ecpu),'-t',Tpath])
    #NOTE that mixer can use a control.RData, which can be provided
    if not no_controls:
        logging.info("Running EXCAVATOR2 DataAnalysis")
        logging.info("Running " + " ".join(['EXCAVATORDataAnalysis.pl','-v','-s', ayaml,'-o',os.path.join(
            os.path.dirname(os.path.abspath(ayaml)), 'output/DataAnalysis_w50k'),'-@',str(ecpu),'-t',Tpath, '-i',os.path.join(
        os.path.dirname(os.path.abspath(pyaml)), 'output/DataPrepare_w50k'), '-e','pooling']))
        subprocess.run(['EXCAVATORDataAnalysis.pl','-v','-s', ayaml,'-o',os.path.join(
            os.path.dirname(os.path.abspath(pyaml)), 'output/DataAnalysis_w50k'),'-@',str(ecpu),'-t',Tpath, '-i',os.path.join(
        os.path.dirname(os.path.abspath(pyaml)), 'output/DataPrepare_w50k'), '-e','pooling'])
    else:
        logging.info("No control samples provided, skipping EXCAVATOR2 DataAnalysis")
    logging.info("EXCAVATOR2 Analysis completed")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run EXCAVATOR2')
    parser.add_argument('-j', '--json', help="Path to the miXer json file", required=True)
    arguments = parser.parse_args()
    with open(arguments.json, 'r') as j:
        config = json.load(j)
    tmp_folder = os.path.join(
        os.path.abspath(config['main_outdir_host']),
        config['exp_id'],
        "_excavator2_output/"
    )
    if os.path.exists(tmp_folder) == False:
       os.makedirs(tmp_folder)
       print("EXCAVATOR2 temp folder created")
    target_path, prepare_yaml_file = create_target_yaml(config, tmp_folder)
    dataPrepare_yaml_file = create_prepare_yaml(config, tmp_folder)
    dataAnalysis_yaml_file = create_analysis_yaml(config, tmp_folder)
    df = pd.read_csv(
        config['sample_list'], sep="\t", dtype=str, skiprows=1, header=None, names=[
        "ID", "bamPath", "Gender", "sampleType"
    ])
    no_controls = False
    if len(df[df['sampleType'].str.lower() == 'c']) == 0:
        no_controls = True
    run_EXCA2(
        tyaml=prepare_yaml_file,
        pyaml=dataPrepare_yaml_file,
        ayaml=dataAnalysis_yaml_file,
        Tpath=target_path,
        ecpu=config['threads'],
        no_controls=no_controls)


    

