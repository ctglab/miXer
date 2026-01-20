#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 22:35:12 2024

@author: elia
"""

import os
import math
import argparse
import re
from datetime import datetime
import logging
import sys
import json

# setup logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler('mixerDataset.log'),
        logging.StreamHandler(sys.stdout)
    ])


def natural_sort_key(chrom):
    return [int(t) if t.isdigit() else t for t in re.split(r'(\d+)', chrom)]

# Order bed entries by chromosome (natural ordering) and then by cnv start (POS)
def vcf_sort_key(e):
    return (natural_sort_key(e['Chr']), int(e['Start']))

#Function to apply Phred scaling
def phred_scale(p_error):
    if p_error > 0:
        return -10 * math.log10(p_error)
    return None

class RecordEntry(dict):
    # that's the basic structure of a record entry that could be inherited for both BED and VCF entries    
    def _write_to_vcf(self, vcf_file):
        vcf_str = f"{self['Chr']}\t{self['pos']}\t{self['identifier']}\t.\t<{self['State']}>\t{self['qual']:.2f}\t{self['filter_value']}\t{self['INFO']}\t{self['format_str']}\t{self['sample_str']}"
        vcf_file.write(str(vcf_str) + "\n")


# Parse command line arguments
parser = argparse.ArgumentParser(description='Convert TSV file to VCF format.')
parser.add_argument('-j', '--json', help="Path to the miXer json file", required=True)
parser.add_argument('-ref', '--reference', type=str, default='unspecified', help='Reference genome version (default: "unspecified")')
parser.add_argument('-hc', '--hc_only', action='store_true', help ='Boolean flag to use only pre-filtered HC>=0.9 TSV windows file for VCF creation. False if not specified.')
args = parser.parse_args()
with open(args.json, 'r') as j:
    config = json.load(j)
# Determine VCF output directory
tsv_dir = os.path.join(
    os.path.abspath(config['main_outdir_host']),
    config['exp_id'],
    "mixer_windows"
)
output_dir = os.path.join(
    os.path.abspath(config['main_outdir_host']),
    config['exp_id'],
    "mixer_vcfs"
) 
# Get all subfolders in the tsv_dir and get sample names
for _, dirs, _ in os.walk(tsv_dir):
    # sanitize directory name to remove _TARGET
    samples = [x.replace("_TARGET", "") for x in dirs]
    break

for i in range(len(dirs)):
    #Get current sample id
    sample = samples[i]
    #Get current sample's TSV file dir
    curr_sample_dir = os.path.join(tsv_dir, dirs[i])
    #Define current sample's VCF output dir
    curr_sample_vcf_output_dir = os.path.join(output_dir, dirs[i])
    
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
                fields_keys = [
                    field.strip() for field in fields if field in [
                        "Chr", "Start", "End", "State", "CN", "ProbCall", "p_error", "Median_NRC", "window_length"]]
                idxs = {field:fields.index(field) for field in fields_keys}
            else:
                values = [fields[idxs[field]] for field in fields_keys]
                bed_entry = RecordEntry()
                for k,v in zip(fields_keys, values):
                    bed_entry[k] = v
                bed_entries.append(bed_entry)
    
    #Natural ordering of bed entries
    bed_entries.sort(key=vcf_sort_key)
    
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
        # now bed entries is a collection of RecordEntry objects
        for i, entry in enumerate(bed_entries[1:], start=1): # skip the header
            entry['pos'] = int(entry['Start']) - 1  # Convert to 0-based index
            entry['ProbCall'] = float(entry['ProbCall'])
            entry['entry_number'] = i
            entry['p_error'] = float(entry['p_error'])
            entry['window_length'] = int(float(entry['window_length']))
            entry['svlen_val'] = entry['window_length']
            entry['cn_value'] = entry['CN'] if entry['CN'] != "4+" else 4  # Convert "4+" to 4
            if entry['ProbCall'] > 0.9:
                entry['filter_value'] = "PASS"
            elif 0.7 <= entry['ProbCall'] <= 0.9:
                entry['filter_value'] = "MediumQual"
            else:
                entry['filter_value'] = "LowQual"
            entry['qual'] = phred_scale(entry['p_error'])
            entry['end_val'] = entry['pos'] + entry['svlen_val']
            entry['INFO'] = f"END={entry['end_val']};IMPRECISE;SVLEN={entry['svlen_val']};SVTYPE=CNV;SVCLAIM=D"
            entry['identifier'] = f"miXer_{entry['State']}_{entry['entry_number']}"
            entry['format_str'] = "GT:CN:MNRC:CS"
            entry['sample_str'] = f".:{entry['cn_value']}:{float(entry['Median_NRC']):.2f}:{float(entry['ProbCall']):.2f}"
            entry._write_to_vcf(vcf)
print("VCF file writing done. MiXer run completed")   