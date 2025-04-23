import pandas as pd
import pysam
import os
import warnings
import argparse
import sys


def check_odd(maleList,femaleList):
    if len(maleList) > len(femaleList):
       print("Found {} M and {} F training samples: creating merged dataset using {} M and F".format(len(maleList),len(femaleList),len(femaleList)))
       nrows = len(femaleList)
       maleList = maleList.iloc[:nrows,:]
    elif len(maleList) < len(femaleList):      
       print("Found {} M and {} F training samples: creating merged dataset using {} M and F".format(len(maleList),len(femaleList),len(maleList)))
       nrows = len(maleList)
       femaleList = femaleList.iloc[:nrows,:] 
    else:
       raise ValueError("Check in the config file that there are both M and F training samples") 
    return maleList,femaleList


def check_divergence(fixed_maleList,fixed_femaleList):
    fixed_maleList = fixed_maleList.rename(columns={'ID':'ID_M','Gender':'Gender_M','MeanCvg':'MeanCvg_M'})[['ID_M','Gender_M','MeanCvg_M']]
    fixed_femaleList = fixed_femaleList.rename(columns={'ID':'ID_F','Gender':'Gender_F','MeanCvg':'MeanCvg_F'})[['ID_F','Gender_F','MeanCvg_F']]
    startLoop = fixed_maleList['ID_M']
    mergedList=pd.DataFrame()        
    for i in range(len(startLoop)):
        for j in range(len(startLoop)):
            match2eval = pd.DataFrame(data={'ID_F': [fixed_femaleList['ID_F'][i]], 'ID_M': [fixed_maleList['ID_M'][j]], 'MeanCvg_F': [fixed_femaleList['MeanCvg_F'][i]], 'MeanCvg_M' : [fixed_maleList['MeanCvg_M'][j]]})
            mergedList = pd.concat([mergedList,match2eval])
    mergedList=mergedList.reset_index(drop=True)
    mergedList['madC']  = mergedList.reset_index(drop=True).mad(axis=1)
    sample2merge=pd.DataFrame() 
    while len(mergedList)!=0:
          list2rm = []
          a = mergedList.loc[mergedList.groupby("ID_F")["madC"].idxmin()].drop_duplicates('ID_M')
          sample2merge = pd.concat([sample2merge,a])
          list2rm = list2rm + a.ID_F.to_list() + a.ID_M.to_list()
          mergedList = mergedList[(~mergedList['ID_F'].isin(list2rm)) & (~mergedList['ID_M'].isin(list2rm))]
    sample2merge = sample2merge.reset_index(drop=True)
    return sample2merge


def mergeNindex(label,final,fP,mP,MCF,th):
     pysam.merge('-@',th,'-f', final, fP, mP)
     pysam.index('-@',th,final)
     updatedConfigFile = pd.DataFrame(data={'ID':[label],'bamPath':[final],'Gender':['MF'],'sampleType':['train']})
     return updatedConfigFile   


#######################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create XXXY samples for training')
    parser.add_argument('-o', '--output_dir', metavar="", help="path to output dir")
    parser.add_argument('-e', '--exp_name', default = "", metavar="", help="experiment or target name")
    parser.add_argument('-c', '--cvg', metavar="", default = None, required = False, help="File with all coverages created with computeAutosomalCvg.py script")
    parser.add_argument('-cf', '--config', metavar="", help="tab separated file with all samples autosomes coverage: ID meancvg")
    parser.add_argument('-th', '--thread', metavar="", help="number of threads to use for merging BAM files", type=str)
    
    if len(sys.argv)==1:
        print()
        print("Usage: python", sys.argv[0]," --help")
        print()
        raise ValueError("No arguments provided")
    arguments = parser.parse_args()
    ########### DECLARE VARIABLES:
    
    ###inputs:
    conf = pd.read_table(arguments.config, names= ['ID','bamPath','Gender','sampleType'], sep="\t")
    conf.sampleType = conf.sampleType.str.replace(' ', '').str.lower()
    conf_df = conf[conf['sampleType'].str.contains(r'train')].reset_index(drop=True)

    ###outputs:
    configBasename = os.path.basename(arguments.config)
    configU = os.path.join(arguments.output_dir,configBasename+'_forExca2')
    conf.to_csv(configU, mode='w', sep='\t', index=False, header=True)

    if not conf_df.empty:
       cvg = pd.read_table(arguments.cvg,names= ['ID','MeanCvg'], sep="\t") 
       ###tmp dir
       absP = os.path.abspath(arguments.output_dir)
       tmp_folder = os.path.join(absP,arguments.exp_name + "_mergedBam_tmp/")
       if os.path.exists(tmp_folder) == False:
          os.makedirs(tmp_folder)
       print("temp folder for MF BAM created")
       ## add mean cvg to config file
       conf_df.Gender = conf_df.Gender.str.replace(' ', '').str.lower()
       merged_list = pd.merge(conf_df,cvg, on='ID')
       ##separate M and F in two dataframes
       bam_M = merged_list.loc[merged_list['Gender']=="m"].reset_index(drop=True)
       bam_F = merged_list.loc[merged_list['Gender']=="f"].reset_index(drop=True)
       if len(bam_M) != len(bam_F):
          bam_M,bam_F = check_odd(bam_M,bam_F) ##reduce dimensionality of the two list of BAM files based on the smallest list   
       sample2merge = check_divergence(bam_M,bam_F)
       outliers = sample2merge.loc[sample2merge['madC']>10]
       if not outliers.empty:
          warnings.warn("The following Males and Female BAM files for training have divergent Mean Coverages: are they from the same sequencing batch? Results can be unpredictable")  
          print('{}\n'.format(outliers))
       
       #start merging samples:   
       start = sample2merge['madC']   
       for i in range(len(start)):
            idMerged = sample2merge['ID_F'][i]+'_'+sample2merge['ID_M'][i]
            outfile = os.path.join(tmp_folder,idMerged +'.bam')  ###crearli in cartella temporanea: servono ancora per exca2
            cvgF = merged_list.loc[merged_list['ID']==sample2merge['ID_F'][i]].MeanCvg.to_list()[0]
            pathF= merged_list.loc[merged_list['ID']==sample2merge['ID_F'][i]].bamPath.to_list()[0]
            pathM= merged_list.loc[merged_list['ID']==sample2merge['ID_M'][i]].bamPath.to_list()[0]
            print('Now merging female sample {} and male sample {}'.format(sample2merge['ID_F'][i],sample2merge['ID_M'][i]))
            updatedConfigFile = mergeNindex(idMerged,outfile,pathF,pathM,cvgF,arguments.thread)
            updatedConfigFile.to_csv(configU, mode='a', sep='\t', index=False, header=False)
    else:
       conf.to_csv(configU, mode='w', sep='\t', index=False, header=True) #save the initial config in the correct format
       print("No training samples found, skipping XXXY creation...") 
       
        



        
