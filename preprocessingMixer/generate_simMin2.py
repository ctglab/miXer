# -*- coding: utf-8 -*-
#simulate class -2 using Rdata of control samples
#this script needs the control samples of a kit, it will simulate chrX XLR exons with a double deletion
#it will do so by setting the NRC data of random XLR exons from the control samples very close to zero (positive values) drawn from a Poisson distribution

#since MeanCvg is needed and as of now is not available for control samples, it will be simulated too
#by estimating its distribution from a pre-made version of the merged XLR training dataset
#the MeanCvg of the control sample will be then randomly drawn from the estimated distribution

#the script will create (n_control_samples)*100 simulated exons for the training if n_control_samples >= 10
#else, it will create 1000 simulated exons

#TODO integrare in merged_maker.py ##### ha senso?

####################### TODO per ora gli serve anche un merged XLR già creato per stimare la distribuzione del suo MeanCvg 
#da quella distribuzione viene simulato il MeanCvg dei campioni di controllo, che per ora non lo hanno, va bene così?
#######################

import sys
from pybedtools import BedTool,helpers,contrib
import pandas as pd
import os
import argparse
import numpy as np
from rpy2.robjects import pandas2ri,r
import rpy2.robjects as robjects
import scipy
import scipy.stats
#import matplotlib.pyplot as plt
from rpy2.robjects.conversion import localconverter

tmp_folder = os.path.join(os.getcwd(),"tmp_folder")
if os.path.exists(tmp_folder) == False:
   os.makedirs(tmp_folder)
helpers.set_tempdir(tmp_folder)

#setting minimum/maximum simulated meancvg values
min_Smcvg = 50
max_Smcvg = 200
# multiplying factor for redacted nrc pool
# based on the increase of cvg, nrc should grow (almost linearly) (as observed when creating the training dataset)
# this should happen even for double deletions
# cvg = min_Smcvg --> factor = 1
# cvg = max_Smcvg --> factor = max_Smcvg/min_Smcvg (e.g. 200/50 = 4)

points = [(min_Smcvg, 1), (max_Smcvg, max_Smcvg / min_Smcvg)]
x_coords, y_coords = zip(*points)
A = np.vstack([x_coords, np.ones(len(x_coords))]).T
m, q = np.linalg.lstsq(A, y_coords, rcond=None)[0]
m = round(m, 4)
q = round(q, 4)

#print("Multiplication Factor function parameters m: {} q: {}".format(m,q))

#TODO resources folder, anche per generate_miXer_datasets.py
#### extract_nrc: function that load and extract a R matrix from RData. Matrix in then converted in a pandas df, and df is formatted to keep only the coordinates and RC of IN-target regions 
def extract_nrc(rdatafile):
    pandas2ri.activate()
    r['load'](rdatafile)
    matrix = robjects.globalenv['MatrixNorm']
    df = r['as.data.frame'](matrix)
    with localconverter(robjects.default_converter + pandas2ri.converter):
         df_conv = robjects.conversion.rpy2py(df)
    df_frt = df_conv.loc[df_conv['Class'] == "IN"].reset_index(drop=True).loc[:, ['chrom','start','end','RCNorm']]
    return df_frt

#### fix_target: keep only chromosomal coordinates (remove additional columns if present)
def fix_target(a):
    target_slice=a.iloc[:,:3]
    if any(c.isalpha() for c in str(target_slice.iloc[0,2])) == True:
       target_slice = target_slice.drop([0]).reset_index(drop=True)
    target_bed = BedTool.from_dataframe(target_slice).sort()
    if not target_bed or target_bed == "0":
       raise ValueError("Check that target file is properly formatted")
    else:
        return target_bed
    
#### annotate_target: add target-specifi features (mappability, GC content, length)--> made for all chromosomes
def annotate_target(tar,ref,mapp_file,gapcen):
    ##GC_content: 
    nuc = BedTool.nucleotide_content(fix_target(tar), fi=ref)
    df_with_gc = pd.read_table(nuc.fn).loc[:, ['#1_usercol', '2_usercol', '3_usercol', '5_pct_gc']]
    df_with_gc_bed = BedTool.from_dataframe(df_with_gc).sort()
    ##Mappability:
    #mapp_file = BedTool(mapp).sort()
    intersect_map = df_with_gc_bed.intersect(mapp_file, wo=True)
    df_with_mapp = pd.read_table(intersect_map.fn, names=['chrom', 'start', 'end', 'gc', 'chr_map', 'start_map', 'end_map', 'mappability', 'count']).loc[:, ['chrom','start','end','gc','mappability']]
    df_with_mapp_bed = BedTool.from_dataframe(df_with_mapp).sort()
    merge_df = df_with_mapp_bed.merge(c=[4,5], o=['distinct','mean'])
    ##remove GAP regions
    final_target=merge_df.intersect(gapcen, v=True)
    final_target_df = pd.read_table(final_target.fn, names=['Chr', 'Start', 'End', 'GC_content', 'Mappability'])
    ##length:
    final_target_df['Length'] = final_target_df['End'].sub(final_target_df['Start'])
    if final_target_df.empty:
       raise ValueError("Check that target file is not empty and properly formatted")
    else:    
       return final_target_df
    
def create_annotated_simXLR(target,fasta,mappability,gapcen,xlrfile):
    
    final_target=BedTool.from_dataframe(annotate_target(target,fasta,mappability,gapcen))
    xlr = BedTool(xlrfile).sort()
    ##file XLR
    xlr_on_target= xlr.intersect(fix_target(target)).sort().merge()
    intersect_xlr = final_target.intersect(xlr_on_target, wo=True)
    df_xlr = pd.read_table(intersect_xlr.fn, header=None).loc[:,:5]
    print(df_xlr.head())
    if df_xlr.empty:
       raise ValueError("Check that xlr file is not empty and properly formatted")
    else:
       print('Creating target file with X-Linked Recessive regions..') 
    return df_xlr

def make_dataset(samplename,infile,pool,cvg,df):
    nrc_file = BedTool.from_dataframe(infile).sort()
    frt_target = BedTool.from_dataframe(df).sort()  
    nrc_pool = BedTool.from_dataframe(pool).sort()
    ##poolF nrc intersection with target
    intersect_pool= nrc_pool.intersect(frt_target)
    df_poolF = pd.read_table(intersect_pool.fn, names=['chrom', 'start', 'end', 'nrc_pool'])
    ##sample nrc annotation with target
    intersect_nrc = nrc_file.intersect(frt_target, wo=True)
    df_with_nrc= pd.read_table(intersect_nrc.fn, names=['chrom', 'start', 'end', 'nrc', 'chr_tar', 'start_tar', 'end_tar', 'gc', 'mappability','length', 'count']).loc[:, ['chrom','start','end','gc','mappability','length','nrc']]
    df_with_nrc['id'] = samplename
    df_with_cvg = pd.merge(df_with_nrc,cvg, on='id')
    df_with_nrcpool= pd.merge(df_with_cvg,df_poolF, on=['chrom','start','end'])
    
    multFactor = df_with_nrcpool["multFactor"].unique()[0]
    #######################qui viene calcolato l'NRC poolnorm, devo cambiare df_with_nrcpool
    #devo sostituire il vettore df_with_nrcpool con uno lungo uguale ma generato a random con una gamma che da valori
    #vicinissimi allo zero
    #plt.figure()
    new_redacted_nrcpool = np.random.beta(1.2, 10, size = df_with_nrcpool.shape[0])*0.01*multFactor
    #plt.hist(new_redacted_nrcpool, density=True, bins=30)
    #plt.show()
    df_with_nrcpool['nrc_poolNorm']= new_redacted_nrcpool/df_with_nrcpool['nrc_pool']
    #######################
    
    #plt.figure()
    df_with_nrcpool['nrc_poolNorm'] = np.log2(df_with_nrcpool['nrc_poolNorm'])
    #df_with_nrcpool["nrc_poolNorm"].hist()
    #plt.show()
    #To simulate class -2 (the only one present now)
    df_with_nrcpool.loc[df_with_nrcpool["gender"]=="MM", "gender"] = -2
    df_final = df_with_nrcpool.loc[:, ['chrom','start','end','gc','mappability','length', 'cvg','nrc_poolNorm', 'id', 'gender']]
    if df_final.empty:
       raise ValueError('Dataset for sample %s is empty' % samplename)
    else: 
       return df_final

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create simulated -2 class for training')
    parser.add_argument('-o', '--output_dir', metavar="", help="path to output dir")
    parser.add_argument('-pmm', '--pre_made_merged', metavar="", help=" path to premade XLR merged dataset for meancvg estimation")
    parser.add_argument('-t', '--target', metavar="", help="target bed file - the same used for alignment and EXCAVATOR2")
    parser.add_argument('-r','--ref', metavar="", help="reference fasta file - the same used for alignment and EXCAVATOR2", type=argparse.FileType('r'))
    parser.add_argument('-p','--par', metavar="", help="Pseudo-autosomal regions file")
    parser.add_argument('-g','--gap', metavar="", help="EXCAVATOR2 GAP regions file")
    parser.add_argument('-cm','--centromeres', metavar="", help="EXCAVATOR2 centromeres file")
    parser.add_argument('-f','--poolF_nrc', metavar="", help="EXCAVATOR2 Rdata file with nrc from pool of females")
    parser.add_argument('-x', '--xlr', metavar="", help="XLR genes file")
    parser.add_argument('-m','--mapp', metavar="", help="EXCAVATOR2 bigwig mappability file")
    ####################### Definisci questo file: idea: carica un merged e stima la distribuzione del MeanCvg
    parser.add_argument('-c', '--coverage_file', metavar="", help="tab separated file with all samples autosomes coverage: ID meancvg", type=argparse.FileType('r'))
    #######################
    parser.add_argument('-e', '--exp_name', metavar="", help="experiment or target name")
    parser.add_argument('-n', '--samples_nrc', nargs='*', metavar="FILES", help="EXCAVATOR2 RCNorm files generated with DataPrepare")
    
    if len(sys.argv)==1:
        print()
        print("Usage: python", sys.argv[0]," --help")
        print()
        raise ValueError("No arguments provided")
    arguments = parser.parse_args()
    
    #load fake coverage file, will change the cvg column with a simulated one
    cvg_df = pd.read_csv(arguments.coverage_file, sep="\t", index_col=False, names=['id', 'cvg'], header = 0) 
    
    #usually, F samples are a subset of the training samples and are found in the same folder
    #only samples specified in cvg_df will be processed
    to_process = cvg_df["id"].unique()
    #load template merged dataset
    model_xlr = pd.read_csv(arguments.pre_made_merged, sep = "\t")
    #fit distribution to MeanCvg
    mu, sigma = scipy.stats.norm.fit(model_xlr.loc[(model_xlr["MeanCvg"] >= min_Smcvg) & (model_xlr["MeanCvg"] <= max_Smcvg)]["MeanCvg"])
    #TODO evaluate limits of meancvg distribution
    #generate n mean coverages for n samples
    fakeMcvgs = []
    
    while len(fakeMcvgs) < len(cvg_df):
        t = np.random.normal(mu, sigma, 1)[0]
        
        if (t >= min_Smcvg) & (t <= max_Smcvg):
            fakeMcvgs.append(t)
    
    cvg_df["gender"] = "MM"
    cvg_df["cvg"] = fakeMcvgs
    cvg_df["multFactor"] = cvg_df["cvg"]*m + q
    
    # ax = plt.subplot(111)
    # ax.hist(model_xlr["MeanCvg"], color='green', density=True)
    # ax.hist(fakeMcvgs, color = 'red', density = True)
    # plt.show()
    
    target_df=pd.read_csv(arguments.target, sep="\t", index_col=False, header=None)
    
    nrc_pool_rdata2bed = extract_nrc(arguments.poolF_nrc)
 
    ###create folders:
    out_dir = arguments.output_dir
    if not os.path.exists(out_dir):
       os.mkdir(out_dir)
    target_dir = os.path.join(out_dir,'targets')
    if not os.path.exists(target_dir):
       os.mkdir(target_dir)
    dataset_dir = os.path.join(out_dir,'datasets')
    if not os.path.exists(dataset_dir):
       os.mkdir(dataset_dir)
    sample_dir = os.path.join(dataset_dir,arguments.exp_name)
    if not os.path.exists(sample_dir):
       os.mkdir(sample_dir)
    
    ###create filenames:
    tarname = os.path.basename(arguments.target)
    xlr_tarOut= os.path.join(target_dir,tarname+'_filtered_min2SimXLR.bed')
    xlr_setOut= os.path.join(sample_dir,'ALL_SAMPLE_min2Sim_XLR.txt.gz')
    
    print('Starting target annotation...')
    checkb37 = target_df[target_df[0].str.startswith(('chr'))]
    ##mappability: from bigwig to bedtool object 
    map_bt = contrib.bigwig.bigwig_to_bedgraph(arguments.mapp) #<- final file
    if checkb37.empty:
       map_df = pd.read_table(map_bt.fn, names=['chr','start','end','mapp'])
       map_df['chr'] = map_df['chr'].str.replace('chr','')
       map_bt = BedTool.from_dataframe(map_df).sort()
      
    ##open centromeres and GAP files and concatenate 
    centro = pd.read_table(arguments.centromeres)
    gap = pd.read_table(arguments.gap)
    #extract only chrom interval
    gap = gap[['chrom','chromStart','chromEnd']]
    centrogap = pd.concat([centro,gap]).sort_values(by=['chrom','chromStart','chromEnd']).drop_duplicates().reset_index(drop=True)
    if checkb37.empty:
       centrogap['chrom'] = centrogap['chrom'].str.replace('chr','')
    centrogap_bt = BedTool.from_dataframe(centrogap).sort() #<- final file
    
    df_xlr = create_annotated_simXLR(target_df,arguments.ref, map_bt,centrogap_bt,arguments.xlr)
    
    df_xlr.to_csv(xlr_tarOut, sep='\t', index=False, header=False)
    outFxlr = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','MeanCvg','NRC_poolNorm', 'ID','Class'])
    outFxlr.to_csv(xlr_setOut, mode='w', sep='\t', index=False)
    
    print('Starting samples annotation...')
    for line in arguments.samples_nrc:
        name = '.'.join(os.path.basename(line).split(".")[:-2])
        if name in to_process:
            sample=line.strip()
            print('Now processing sample',name)
            filename = os.path.basename(sample).split(".")[0]
            nrc_rdata2bed = extract_nrc(sample)    
            make_dataset(filename,nrc_rdata2bed,nrc_pool_rdata2bed,cvg_df,df_xlr).to_csv(xlr_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            #non ho bisogno di nsd e sd, solo degli XLR per il training visto che li sto simulando per il training

##note:
###dare una output folder dove poter mettere tutti i dataset (dataset --> target_name/Experiment_name --> file)

