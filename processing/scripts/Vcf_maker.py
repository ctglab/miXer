#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:35:12 2024

@author: elia
"""

import os
import math
import argparse
from datetime import datetime

#Function to apply Phred scaling
def phred_scale(p_error):
    if p_error > 0:
        return -10 * math.log10(p_error)
    return None

# Parse command line arguments
parser = argparse.ArgumentParser(description='Convert TSV file to VCF format.')
parser.add_argument('-td','--tsv_dir', type=str, required=True, help='Directory containing TSV files from HMM output')
parser.add_argument('-od', '--output_dir', type=str, help='Output directory for VCF files', default = None)
parser.add_argument('-ref', '--reference', type=str, default='unspecified', help='Reference genome version (default: "unspecified")')
parser.add_argument('-hc', '--hc_only', action='store_true', help ='Boolean flag to use only pre-filtered HC>=0.9 TSV windows file for VCF creation. False if not specified.')
args = parser.parse_args()

# Determine VCF output directory

if args.output_dir is None:
    bed_parent_dir = os.path.dirname(os.path.dirname(args.tsv_dir)) #TODO fix it?
    vcf_dir = os.path.join(bed_parent_dir, "VCF")
else:
    vcf_dir = args.output_dir

# Get all subfolders in the tsv_dir and get sample names
for _, dirs, _ in os.walk(args.tsv_dir):
    # sanitize directory name to remove _TARGET
    samples = [x.replace("_TARGET", "") for x in dirs]
    break

for i in range(len(dirs)):
    #Get current sample id
    sample = samples[i]
    #Get current sample's TSV file dir
    curr_sample_dir = os.path.join(args.tsv_dir, dirs[i])
    #Define current sample's VCF output dir
    curr_sample_vcf_output_dir = os.path.join(vcf_dir, dirs[i])
    
    #Check if output dir exists, if not, create it
    if not os.path.isdir(curr_sample_vcf_output_dir):
        os.makedirs(curr_sample_vcf_output_dir)
    #Get correct TSV file
    #Get all files in directory
    file_list = os.listdir(curr_sample_dir)
    
    if args.hc_only:
        #Use file with only high confidence windows
        bed_file = [filename for filename in file_list if 'PASS' in filename][0]
        
    else:
        #Use file with all windows
        bed_file = [filename for filename in file_list if 'PASS' not in filename][0]
    
    
    # Read the TSV file
    bed_entries = []
    with open(os.path.join(curr_sample_dir, bed_file), 'r') as f:
        for line_num, line in enumerate(f):
            fields = line.strip().split('\t')
            if line_num == 0:
                # The first line is the header; find the index of Median_NRC
                chr_idx = fields.index("Chr")
                start_idx = fields.index("Start")
                end_idx = fields.index("End")
                state_idx = fields.index("State")
                cn_idx = fields.index("CN")
                probcall_idx = fields.index("ProbCall")
                p_error_idx = fields.index("p_error")
                median_nrc_idx = fields.index("Median_NRC")
                window_length_idx = fields.index("window_length")
                continue  # skip header
            bed_entries.append(fields)
        
    #Get current date
    current_date = datetime.now().strftime('%Y%m%d')
    
    #Create VCF header
    header_lines = [
        "##fileformat=VCFv4.4",
        f"##fileDate={current_date}",
        "##source=miXer_v1.0",
        f"##reference={args.reference}",
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##FILTER=<ID=PASS,Description="High-quality calls with CS higher than 0.9">',
        '##FILTER=<ID=MediumQual,Description="Calls with CS between 0.7 and 0.9">',
        '##FILTER=<ID=LowQual,Description="Calls with CS lower than 0.7">',
        '##QUAL=<ID=PQS,Number=1,Type=Float,Description="Phred-scaled quality score of the probability of a wrong CNV type assignment">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant described in this record">',
        '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
        '##INFO=<ID=SVCLAIM,Number=1,Type=String,Description="D for claiming read depth based structural variant calling">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Window Genotype">',
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number inferred by majority voting of SVM calls over encompassing target regions">',
        '##FORMAT=<ID=MNRC,Number=1,Type=Float,Description="Median NRC_PoolNorm over encompassing target regions">',
        '##FORMAT=<ID=CS,Number=1,Type=Float,Description="Confidence score of Window copy-number call">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    ]

    #Join header lines with newline characters
    header = "\n".join(header_lines) + "\n"
    
    #write to VCF file
    curr_sample_vcf_output = os.path.join(curr_sample_vcf_output_dir, sample) + ".vcf"
    
    with open(curr_sample_vcf_output, 'w') as vcf:
        vcf.write(header)
        
        for entry_number, entry in enumerate(bed_entries, start=1):
            #CHROM entry
            chrom = entry[chr_idx]
            #POS entry <- Position of the basis before the event
            pos = int(entry[start_idx]) - 1
            #End 
            end = entry[end_idx]
            #State <- DEL/DUP, corresponds to HMM state #TODO check for DDEL -> DEL
            state = entry[state_idx]
            #Call probability
            probcall = float(entry[probcall_idx])
            #Non-phred scaled error probability
            p_error = float(entry[p_error_idx])
            # Window length
            window_length = int(entry[window_length_idx])
            # SVLEN value
            svlen_val = window_length
            # Copy number (CN) value
            cn_value = entry[cn_idx]
            #Check if copy number is 4+ (meaning 4 or more copies of a TR are present, as of now we cannot discriminate between 4 and more than 4 copies)
            #If it is 4+ set to 4
            if cn_value == "4+":
                cn_value = 4
            # NRC_PoolNorm value
            nrc_poolNorm = float(entry[median_nrc_idx])
            
            #Set FILTER based on ProbCall
            if probcall > 0.9:
                filter_value = "PASS"
            elif 0.7 <= probcall <= 0.9:
                filter_value = "MediumQual"
            else:
                filter_value = "LowQual"
            
            # Apply Phred scaling
            qual = phred_scale(p_error)
            
            # Calculate END value for INFO field
            end_val = pos + svlen_val
            
            # Prepare INFO field
            info = f"END={end_val};IMPRECISE;SVLEN={svlen_val};SVTYPE=CNV;SVCLAIM=D"
            
            # ID column value
            identifier = f"miXer_{state}_{entry_number}"
            
            # Prepare VCF line
            #REF is hardcoded to "." since it should be the basis (ACTG) before the event # TODO add in next release?
            format_str = "GT:CN:MNRC:CS"
            sample_str = f".:{cn_value}:{nrc_poolNorm:.2f}:{probcall:.2f}"
            
            vcf_line = f"{chrom}\t{pos}\t{identifier}\t.\t<{state}>\t{qual:.2f}\t{filter_value}\t{info}\t{format_str}\t{sample_str}\n"

            vcf.write(vcf_line)

print("VCF file writing done. MiXer run completed")   