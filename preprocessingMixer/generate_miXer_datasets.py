import sys
from pybedtools import BedTool,helpers,contrib
import pandas as pd
import os
import argparse
import numpy as np
from rpy2.robjects import pandas2ri,r
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
from datetime import datetime
import glob
import shutil
import json
import gc
import logging

# setup logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler('mixerDataset.log'),
        logging.StreamHandler(sys.stdout)
    ])


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
def fix_target(a: pd.DataFrame) -> BedTool:
    target_slice=a.iloc[:,:3]
    if any(c.isalpha() for c in str(target_slice.iloc[0,2])) == True:
       target_slice = target_slice.drop([0]).reset_index(drop=True)
    target_bed = BedTool.from_dataframe(target_slice).sort()
    if not target_bed or target_bed == "0":
       raise ValueError("Check that target file is properly formatted")
    else:
        return target_bed
        
#### annotate_target: add target-specifi features (mappability, GC content, length)--> made for all chromosomes
def annotate_target(tar: pd.DataFrame, ref: str, mapp: str, gapcen: BedTool):
    target_bed = fix_target(tar)
    # pybedtools could do things in pipe. 
    final_bed = (
        # gc content
        target_bed.nucleotide_content(fi=ref)
        # pybedtools uses 0-based indexing for cut, so we use 0, 1, 2, 4.
        .cut([0, 1, 2, 4])
        # map is more efficient that intersect - pandas merge
        .map(b=mapp, c=4, o='mean')
        # remove gaps
        .intersect(gapcen, v=True)
    )
    # convert then to pandas
    final_target = final_bed.to_dataframe(
        names=['Chr', 'Start', 'End', 'GC_content', 'Mappability']
    )
    final_target['Length'] = final_target['End'].sub(final_target['Start'])
    logging.info(f"Annotation completed")
    # Check for an empty result after all filtering.
    if final_target.empty:
       raise ValueError("Check that target file is not empty and properly formatted after filtering.")
    else:    
       return final_target

#### create_training_datasets function: select only chrX and split the chromosome in the XLR, SD and NSD datasets
def create_training_datasets(targetWfeatures,target,pseud,xlrfile,segdupfile):    
    ##Select only chrX and convert df in pybedtool object
    df_with_chrX=targetWfeatures.loc[(targetWfeatures.chrom == "chrX") | (targetWfeatures.chrom == "X") | (targetWfeatures.chrom == "ChrX")]
    df_with_chrX_bt = BedTool.from_dataframe(df_with_chrX).sort()
    ##remove Pseudo-Autosomal Regions
    final_target=df_with_chrX_bt.intersect(pseud, v=True)
    ##select XLR target regions
    xlr_on_target= xlrfile.intersect(fix_target(target)).sort().merge()
    intersect_xlr = final_target.intersect(xlr_on_target, wo=True)
    #df_xlr = pd.read_table(intersect_xlr.fn, header=None).loc[:,:5]
    df_xlr = pd.read_table(intersect_xlr.fn, names=['chrom', 'start', 'end', 'gc', 'mappability', 'length', 'nrc_poolnorm', 'id', 'countX', 'chr_seg', 'start_Seg', 'end_seg', 'count']).loc[:, ['chrom','start','end','gc','mappability','length','nrc_poolnorm','id']].drop_duplicates()
    if df_xlr.empty:
       raise ValueError("Check that xlr file is not empty and properly formatted")
    else:
       logging.info('Creating target file with X-Linked Recessive regions..') 
    intersect_no_xlr = final_target.intersect(xlr_on_target, v=True)
    ##select no_SegDup target regions
    intersect_no_segdup= intersect_no_xlr.intersect(segdupfile, v=True)
    df_nosegdup =pd.read_table(intersect_no_segdup.fn, header=None)
    if df_nosegdup.empty:
       raise ValueError("Check that segdup file is not empty and properly formatted")
    else: 
       logging.info('Creating target file with Normal regions..')
    ##select SegDup target regions
    intersect_segdup= intersect_no_xlr.intersect(segdupfile, wo=True)
    df_with_segdup = pd.read_table(intersect_segdup.fn, names=['chrom', 'start', 'end', 'gc', 'mappability', 'length', 'nrc_poolnorm', 'id', 'countX', 'chr_seg', 'start_Seg', 'end_seg', 'count']).loc[:, ['chrom','start','end','gc','mappability','length','nrc_poolnorm','id']].drop_duplicates()
    if df_nosegdup.empty:
       raise ValueError("Check that segdup file is not empty and properly formatted")
    else: 
       logging.info('Creating target file with Segmental Duplication regions..') 
    return df_xlr, df_nosegdup, df_with_segdup

#### make_dataset function: add sample-specific features to the dataset (used for both training and calling)
def make_dataset(samplename,infile,pool,df,is_mf =False, mf_mapp_threshold = 0.99):
    nrc_file = BedTool.from_dataframe(infile).sort()
    frt_target = BedTool.from_dataframe(df).sort()  
    nrc_pool = BedTool.from_dataframe(pool).sort()
    ##poolF nrc intersection with target
    intersect_pool= nrc_pool.intersect(frt_target)
    df_poolF = pd.read_table(intersect_pool.fn, names=['chrom', 'start', 'end', 'nrc_pool'])
    ##sample nrc annotation with target
    intersect_nrc = nrc_file.intersect(frt_target, wo=True)
    df_with_nrc= pd.read_table(intersect_nrc.fn, names=['chrom', 'start', 'end', 'nrc', 'chr_tar', 'start_tar', 'end_tar', 'gc', 'mappability','length', 'count']).loc[:, ['chrom','start','end','gc','mappability','length', 'nrc']]
    df_with_nrc['id'] = samplename
    #df_with_cvg = pd.merge(df_with_nrc,cvg, on='id')
    df_with_nrcpool= pd.merge(df_with_nrc,df_poolF, on=['chrom','start','end'])
    df_with_nrcpool['nrc_poolNorm']= df_with_nrcpool['nrc']/df_with_nrcpool['nrc_pool']
    df_with_nrcpool['nrc_poolNorm'] = np.log2(df_with_nrcpool['nrc_poolNorm'])

    #if not is_mf:
    #    nrc_median = df_with_nrcpool[~df_with_nrcpool["chrom"].isin(["chrX", "ChrX", "chrx", "X", "chrY", "ChrY", "chry", "Y"])]["nrc_poolNorm"].median()
    
    nrc_median = df_with_nrcpool[~df_with_nrcpool["chrom"].isin(["chrX", "ChrX", "chrx", "X", "chrY", "ChrY", "chry", "Y"])]["nrc_poolNorm"].median()
    
    #df_with_nrcpool['nrc_poolNorm'] = df_with_nrcpool['nrc_poolNorm'] - nrc_median
        
    df_final = df_with_nrcpool.loc[:, ['chrom','start','end','gc','mappability','length', 'nrc_poolNorm','id']]
    if df_final.empty:
       raise ValueError('Dataset for sample %s is empty, check first config file' % samplename)
    else: 
       if is_mf:
           logging.info("Selecting High mappability (>= {}) Autosomal exons to train the SVM to detect double (or more) duplications".format(mf_mapp_threshold))
           no_chrx = df_final.loc[~df_final["chrom"].isin(["chrX", "ChrX", "X", "x"])]
           #compatibility with make training dataset columns
           no_chrx.columns = ['chrom','start','end','gc','mappability','length','nrc_poolnorm', 'id']
           no_chrx_high_mapp = no_chrx.loc[no_chrx["mappability"] > mf_mapp_threshold]
           no_chrx_high_mapp_median = no_chrx_high_mapp["nrc_poolnorm"].median()
           no_chrx_high_mapp = no_chrx_high_mapp.loc[no_chrx_high_mapp["nrc_poolnorm"] > no_chrx_high_mapp_median]
           # print(no_chrx.shape)
           # print(no_chrx_high_mapp.shape)
           #print(no_chrx_high_mapp["nrc_poolnorm"].median())
           
           return df_final, nrc_median, no_chrx_high_mapp
       else:
           return df_final, nrc_median

#####add gender based on config file
def addclass(training_df, gen, filename):   
    if gen.lower() == "m":
       training_df['Class'] = "-1"
    elif gen.lower() == "f":
         training_df['Class'] = "0"  
    elif (gen.lower() == "mf") | (gen.lower() == "fm"):
         training_df['Class'] = "1"
    else:
        raise ValueError('check gender in config file for sample' % filename)
    return training_df 

def addclass_ddup(ddup_df):
    ddup_df["Class"] = "2"
    return ddup_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create annotated datasets for training')
    parser.add_argument('-j', '--json', help="Path to the miXer json file", required=True)
    parser.add_argument('-s', '--samples', metavar="", required=True, help="miXer samples sheet")
    parser.add_argument('-x', '--xlr', metavar="", help="OPTIONAL: XLR genes file", required=False)
    parser.add_argument('-sd', '--segdup', metavar="", help="OPTIONAL: Segmental Duplications regions file", required=False)
    arguments = parser.parse_args()
    with open(arguments.json, 'r') as j:
        config = json.load(j)
    out_dir = os.path.join(
        os.path.abspath(config['main_outdir_host']),
        config['exp_id']
    )
    # While this is so hardcoded, these are default paths created during preprocessing
    poolF_nrc = os.path.join(
        os.path.abspath(config['main_outdir_host']),
        config['exp_id'],
        "_excavator2_output",
        "output",
        "DataAnalysis_w50k",
        "Control",
        "RCNorm",
        "Control.NRC.RData"
    )
    samples_nrc = glob.glob(os.path.join(
        os.path.abspath(config['main_outdir_host']),
        config['exp_id'],
        "_excavator2_output",
        "output",
        "DataPrepare_w50k",
        "*",
        "RCNorm",
        "*RData"
    ))
    if len(sys.argv)==1:
        print()
        print("Usage: python", sys.argv[0]," --help")
        print()
        raise ValueError("No arguments provided")
    arguments = parser.parse_args()

    ########### DECLARE VARIABLES:
    ####inputs:
    expname = config['exp_id']
    nrc_pool_females = extract_nrc(poolF_nrc)
    #cvg_df = pd.read_csv(arguments.coverage_file, sep="\t", index_col=False, names=['id', 'cvg'])
    #cvg_df.id = cvg_df.id.astype(str)
    target_df = pd.read_csv(config['target'], sep="\t", index_col=False, header=None)
    #Ensuring that first column of target df is a string
    target_df[0] = target_df[0].astype(str)
    checkb37 = target_df[target_df[0].str.startswith(('chr'))]
    caseNtrain_df = pd.read_table(arguments.samples, sep="\t")
    caseNtrain_df.sampleType = caseNtrain_df.sampleType.str.replace(' ', '').str.lower()
    caseNtrain_df.ID = caseNtrain_df.ID.astype(str)
    ##select from config files sample to process for both training and calling
    both = caseNtrain_df[(caseNtrain_df['sampleType'].str.contains(r't')) & (caseNtrain_df['sampleType'].str.contains(r'train'))]
    ## training only
    train_only= caseNtrain_df[caseNtrain_df['sampleType'].isin(['train'])]
    ## calling only
    call_only = caseNtrain_df[caseNtrain_df['sampleType'].isin(['t'])]
    #separo per il training gli M,F only dagli MF
    males_females = train_only.loc[train_only["Gender"].isin(["M", "m", "male", "F", "f", "female", "X", "x"])][["ID", "Gender"]]
    #tengo una lista con gli ID delle femmine usate per fare gli MF simulati
    mfs = train_only.loc[~train_only["Gender"].isin(["M", "m", "male", "F", "f", "female", "X", "x"])][["ID", "Gender"]]
    
    #full list of M, F samples that are used for training = train_only + both
    all_train = pd.concat([train_only, both])
    
    simF = []
    for item in mfs["ID"].to_list():
        found_fem = False
        for thing in males_females["ID"].to_list():
            if thing in item:
                try:
                    gender = all_train.loc[train_only['ID']== thing].Gender.to_list()[0]
                    if gender.lower() == "f":
                        simF.append(thing)
                        found_fem = True
                except:
                    logging.error("Gender not found for {} sample. Check config and coverage file.".format(thing))
                    exit()
        if found_fem is False:
            logging.error("{}'s F sample not found. Check config file".format(item))
            exit()
    
    mfs["simF_ID"] = simF
    ####create folders:
    date = datetime.now().strftime("%Y%m%d-%H%M%S")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    #create pybedtools tmp dir 
    tmp_folder = os.path.join(out_dir,expname+"_miXer_tmp/")
    if os.path.exists(tmp_folder) == False:
        os.makedirs(tmp_folder)
    helpers.set_tempdir(tmp_folder)
    ####EXCAVATOR2 support files 
    ##mappability: from bigwig to bedtool object 
    map_bt = contrib.bigwig.bigwig_to_bedgraph(config['map']) #<- final file
    if checkb37.empty:
       map_df = pd.read_table(map_bt.fn, names=['chr','start','end','mapp'])
       map_df['chr'] = map_df['chr'].str.replace('chr','')
       map_bt = BedTool.from_dataframe(map_df).sort()
    ##open centromeres and GAP files and concatenate 
    centro = pd.read_table(config['centro'])
    gap = pd.read_table(config['gap'])
    #extract only chrom interval
    gap = gap[['chrom','chromStart','chromEnd']]
    centrogap = pd.concat([centro,gap]).sort_values(by=['chrom','chromStart','chromEnd']).drop_duplicates().reset_index(drop=True)
    if checkb37.empty:
       centrogap['chrom'] = centrogap['chrom'].str.replace('chr','')
    centrogap_bt = BedTool.from_dataframe(centrogap).sort() #<- final file
    ##EXCAVATOR2 genome
    genomeFasta = config['ref']
    
    ####miXer-only support files:
    if ((not train_only.empty) | (not both.empty)):
       par_bt = BedTool(config['par']).sort()
       xlr_bt = BedTool(arguments.xlr).sort()
       segdup_bt = BedTool(arguments.segdup).merge().sort()
    
    ####create folders:
    #create training/test set directory
    if ((not train_only.empty) | (not both.empty)):
        dataset_dir = os.path.join(out_dir, f'_datasets_{date}')
        if not os.path.exists(dataset_dir):
           os.makedirs(dataset_dir)
    #create calling dataset directory
    if ((not call_only.empty) | (not both.empty)):
       dataset_test_dir = os.path.join(out_dir, f'_datasets_testing_{date}')
       if not os.path.exists(dataset_test_dir):
          os.makedirs(dataset_test_dir)
   
    ###create filenames and empty dataframes:
    if ((not train_only.empty) | (not both.empty)):
       xlr_setOut= os.path.join(dataset_dir,'ALL_SAMPLE_XLR.txt.gz')
       seg_setOut= os.path.join(dataset_dir,'ALL_SAMPLE_SegDup.txt.gz')     
       noseg_setOut= os.path.join(dataset_dir,'ALL_SAMPLE_noSegDup.txt.gz')
       ddup_setOut=os.path.join(dataset_dir, "ALL_SAMPLE_DDUP.txt.gz")
       outFxlr = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID','Class'])
       outFxlr.to_csv(xlr_setOut, mode='w', sep='\t', index=False)
       outFseg = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID','Class'])
       outFseg.to_csv(seg_setOut, mode='w', sep='\t', index=False)
       outFnoseg = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID','Class'])
       outFnoseg.to_csv(noseg_setOut, mode='w', sep='\t', index=False) 
       outFddup = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID','Class'])
       outFddup.to_csv(ddup_setOut, mode='w', sep='\t', index=False)
       outF_xlr_ddup = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID','Class'])
       
    samples_paths = samples_nrc
    
    #separing MF file paths and other samples file paths
    mf_only = [x for x in samples_paths if '.'.join(os.path.basename(x).split(".")[:-2]) in mfs["ID"].to_list()]
    the_others = [x for x in samples_paths if '.'.join(os.path.basename(x).split(".")[:-2]) not in mfs["ID"].to_list()]
    
    
    ############DATASETS CREATION
    
    logging.info('Starting target annotation...')
    final_target_df=annotate_target(target_df,genomeFasta,map_bt,centrogap_bt)
    ##if there are any training samples, do chrX annotation
    
    logging.info('Starting M and F datasets creation...')
    nrc_median_F_df = pd.DataFrame(columns = ["F_ID", "AutNRC_poolNorm_median"])
    for file in the_others:
        name = '.'.join(os.path.basename(file).split(".")[:-2])
        
        if (name in train_only.ID.astype(str).to_list())  or (name in call_only.ID.astype(str).to_list()) or (name in both.ID.astype(str).to_list()):
            logging.info("Now processing sample {}".format(name))
            nrc_rdata2bed = extract_nrc(file)
            
            full_annot_target, nrc_median = make_dataset(name,nrc_rdata2bed,nrc_pool_females,final_target_df)

            ##samples only required for training
            if name in train_only.ID.astype(str).to_list():
                df_xlr, df_nosegdup, df_with_segdup = create_training_datasets(full_annot_target,target_df,par_bt,xlr_bt,segdup_bt) 
                logging.info('Saving training sample {}\n'.format(name))
                gender = train_only.loc[train_only['ID']== name].Gender.to_list()[0]
                if gender in ["F", "f", "female"]:
                    nrc_median_F_df = pd.concat([nrc_median_F_df, pd.DataFrame(data = { "F_ID": [name], "AutNRC_poolNorm_median": [nrc_median]})])
                    
                addclass(df_xlr,gender,name).to_csv(xlr_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
                addclass(df_with_segdup,gender,name).to_csv(seg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
                addclass(df_nosegdup,gender,name).to_csv(noseg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            ##samples only required for calling   
            elif name in call_only.ID.astype(str).to_list():

                logging.info('Saving use-case processed sample {}\n'.format(name))
                logging.info("Saving sample to {}\n".format(dataset_test_dir))
                target_setOut= os.path.join(dataset_test_dir,name+'_TARGET.txt.gz') 
                outFtarget = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm','ID'])
                outFtarget.to_csv(target_setOut, mode='w', sep='\t', index=False)
                full_annot_target.to_csv(target_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            ##samples to train and call
            elif name in both.ID.astype(str).to_list():
                #no median normalization for training data
                df_xlr, df_nosegdup, df_with_segdup = create_training_datasets(full_annot_target,target_df,par_bt,xlr_bt,segdup_bt) 
                print('Saving training and use-case sample {}\n'.format(name))
                gender = both.loc[both['ID']== name].Gender.to_list()[0]
                if gender in ["F", "f", "female"]:
                    nrc_median_F_df = pd.concat([nrc_median_F_df, pd.DataFrame(data = { "F_ID": [name], "AutNRC_poolNorm_median": [nrc_median]})])
                 
                #Training files
                addclass(df_xlr,gender,name).to_csv(xlr_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
                addclass(df_with_segdup,gender,name).to_csv(seg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
                addclass(df_nosegdup,gender,name).to_csv(noseg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
                #Usecase files
                target_setOut= os.path.join(dataset_test_dir,name+'_TARGET.txt.gz') 
                outFtarget = pd.DataFrame(columns = ['Chr','Start','End','GC_content','Mappability','Length','NRC_poolNorm', 'ID'])
                outFtarget.to_csv(target_setOut, mode='w', sep='\t', index=False)
                full_annot_target.to_csv(target_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
        else:
            continue 
    logging.info('Starting MF datasets creation...')
    for file in mf_only:
        name = '.'.join(os.path.basename(file).split(".")[:-2])
        
        if name in train_only.ID.astype(str).to_list():
            logging.info("Now processing sample {}".format(name))
            #match with f sample tag
            f_samp = mfs.loc[mfs["ID"] == name]["simF_ID"].values[0]
            nrc_rdata2bed = extract_nrc(file)
            full_annot_target, nrc_median, df_ddup = make_dataset(name,nrc_rdata2bed,nrc_pool_females,final_target_df, is_mf = True)
            
            df_xlr, df_nosegdup, df_with_segdup = create_training_datasets(full_annot_target,target_df,par_bt,xlr_bt,segdup_bt)
            
            if df_ddup.shape[0] > int(df_xlr.shape[0]/3): #numero arbitrario, modificabile
                df_ddup = df_ddup.sample(n = int(df_xlr.shape[0]/3), random_state = 42, ignore_index = True)
                
            logging.info('Saving training sample {}\n'.format(name))
            gender = train_only.loc[train_only['ID']== name].Gender.to_list()[0]
            
            #add classes
            df_xlr = addclass(df_xlr,gender,name)
            df_with_segdup = addclass(df_with_segdup,gender,name)
            df_nosegdup = addclass(df_nosegdup,gender,name)
            df_ddup = addclass_ddup(df_ddup)
            
            df_xlr.to_csv(xlr_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            df_with_segdup.to_csv(seg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            df_nosegdup.to_csv(noseg_setOut, mode='a', sep='\t', index=False, header=False, compression="gzip")
            df_ddup.to_csv(ddup_setOut, mode='a', sep='\t', index=False, header=False, compression='gzip')

####clean up
gc.collect()
####rm tmp dir
shutil.rmtree(tmp_folder, ignore_errors=True) 

