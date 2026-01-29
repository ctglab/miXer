#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate inference datasets for miXer SVM prediction.
Creates annotated target region files for each test sample.

This is a refactored version of generate_miXer_datasets.py that removes
all training-related functionality and focuses solely on inference dataset
generation.
"""

import sys
import os
import argparse
import glob
import json
import gc
import shutil
import logging
from datetime import datetime
from typing import Optional

import numpy as np
import pandas as pd
from pybedtools import BedTool, helpers, contrib
from rpy2.robjects import pandas2ri, r
import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter

# --- Constants ---
REQUIRED_CONFIG_KEYS = ['main_outdir_host', 'exp_id', 'target', 'ref', 'map', 'centro', 'gap']
OUTPUT_COLUMNS = ['Chr', 'Start', 'End', 'GC_content', 'Mappability', 'Length', 'NRC_poolNorm', 'ID']
EPSILON = 1e-10  # For log2 stability to avoid division by zero

# Module-level pandas2ri activation (do this once, not per function call)
pandas2ri.activate()

# Setup logger (will be reconfigured in main() with proper output path)
logger = logging.getLogger(__name__)


def setup_logging(log_dir: Optional[str] = None) -> None:
    """
    Configure logging with both file and console handlers.

    Args:
        log_dir: Directory for log file. If None, logs to current directory.
    """
    log_path = os.path.join(log_dir, 'mixerDataset.log') if log_dir else 'mixerDataset.log'

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ],
        force=True  # Override any existing configuration
    )


def validate_config(config: dict) -> None:
    """
    Validate that all required config keys are present.

    Args:
        config: Configuration dictionary loaded from JSON

    Raises:
        ValueError: If required keys are missing
    """
    missing_keys = [key for key in REQUIRED_CONFIG_KEYS if key not in config]
    if missing_keys:
        raise ValueError(f"Missing required config keys: {missing_keys}")

    # Validate that required files exist
    file_keys = ['target', 'ref', 'map', 'centro', 'gap']
    for key in file_keys:
        path = config[key]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Config key '{key}' points to non-existent file: {path}")


def extract_nrc(rdatafile: str) -> pd.DataFrame:
    """
    Load and extract NRC matrix from EXCAVATOR2 RData file.

    Args:
        rdatafile: Path to RData file

    Returns:
        DataFrame with columns: chrom, start, end, RCNorm
        for IN-target regions only.
    """
    r['load'](rdatafile)
    matrix = robjects.globalenv['MatrixNorm']
    df = r['as.data.frame'](matrix)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        df_conv = robjects.conversion.rpy2py(df)
    df_frt = df_conv.loc[df_conv['Class'] == "IN"].reset_index(drop=True)
    return df_frt.loc[:, ['chrom', 'start', 'end', 'RCNorm']]


def fix_target(target_df: pd.DataFrame) -> BedTool:
    """
    Keep only chromosomal coordinates (chr, start, end) from target.

    Args:
        target_df: Target regions DataFrame

    Returns:
        Sorted BedTool object.

    Raises:
        ValueError: If target file is improperly formatted
    """
    target_slice = target_df.iloc[:, :3]
    # Check if first row has alphabetic chars in third column (header row)
    if any(c.isalpha() for c in str(target_slice.iloc[0, 2])):
        target_slice = target_slice.drop([0]).reset_index(drop=True)
    target_bed = BedTool.from_dataframe(target_slice).sort()
    if not target_bed or target_bed == "0":
        raise ValueError("Check that target file is properly formatted")
    return target_bed


def annotate_target(target_df: pd.DataFrame, ref: str, mapp: BedTool,
                    gapcen: BedTool) -> pd.DataFrame:
    """
    Add target-specific features: GC content, mappability, and length.

    Args:
        target_df: Target regions DataFrame
        ref: Path to reference FASTA
        mapp: Mappability BedTool
        gapcen: Gap and centromere regions BedTool

    Returns:
        Annotated target DataFrame with columns:
        Chr, Start, End, GC_content, Mappability, Length
    """
    target_bed = fix_target(target_df)
    final_bed = (
        target_bed.nucleotide_content(fi=ref)
        .cut([0, 1, 2, 4])  # pybedtools uses 0-based indexing
        .map(b=mapp, c=4, o='mean')
        .intersect(gapcen, v=True)  # remove gaps
    )
    final_target = final_bed.to_dataframe(
        names=['Chr', 'Start', 'End', 'GC_content', 'Mappability']
    )
    final_target['Length'] = final_target['End'].sub(final_target['Start'])
    logging.info("Target annotation completed")

    if final_target.empty:
        raise ValueError("Target file is empty after filtering. Check input files.")
    return final_target


def make_dataset(samplename: str, nrc_df: pd.DataFrame,
                 pool_df: pd.DataFrame, target_df: pd.DataFrame,
                 pool_bed: Optional[BedTool] = None,
                 target_bed: Optional[BedTool] = None) -> pd.DataFrame:
    """
    Create sample-specific annotated dataset with NRC_poolNorm.

    Args:
        samplename: Sample identifier
        nrc_df: Sample NRC data from RData
        pool_df: Pool (control) NRC data
        target_df: Annotated target regions
        pool_bed: Pre-computed BedTool from pool_df (optional, for performance)
        target_bed: Pre-computed BedTool from target_df (optional, for performance)

    Returns:
        DataFrame with columns: Chr, Start, End, GC_content, Mappability, Length, NRC_poolNorm, ID
    """
    nrc_bed = BedTool.from_dataframe(nrc_df).sort()

    # Use pre-computed BedTool objects if provided, otherwise create them
    if target_bed is None:
        target_bed = BedTool.from_dataframe(target_df).sort()
    if pool_bed is None:
        pool_bed = BedTool.from_dataframe(pool_df).sort()

    # Pool NRC intersection with target
    intersect_pool = pool_bed.intersect(target_bed)
    df_pool = pd.read_table(intersect_pool.fn, names=['chrom', 'start', 'end', 'nrc_pool'])

    # Sample NRC annotation with target
    intersect_nrc = nrc_bed.intersect(target_bed, wo=True)
    df_with_nrc = pd.read_table(
        intersect_nrc.fn,
        names=['chrom', 'start', 'end', 'nrc', 'chr_tar', 'start_tar', 'end_tar',
               'gc', 'mappability', 'length', 'count']
    ).loc[:, ['chrom', 'start', 'end', 'gc', 'mappability', 'length', 'nrc']]

    df_with_nrc['id'] = samplename

    # Merge with pool and compute normalized NRC
    df_merged = pd.merge(df_with_nrc, df_pool, on=['chrom', 'start', 'end'])

    # Handle division by zero: add small epsilon to avoid inf values
    df_merged['NRC_poolNorm'] = np.log2(
        (df_merged['nrc'] + EPSILON) / (df_merged['nrc_pool'] + EPSILON)
    )
    # Replace any remaining inf values with NaN
    df_merged['NRC_poolNorm'] = df_merged['NRC_poolNorm'].replace([np.inf, -np.inf], np.nan)

    # Rename columns to match expected output format (downstream expects capitalized names)
    df_final = df_merged.rename(columns={
        'chrom': 'Chr',
        'start': 'Start',
        'end': 'End',
        'gc': 'GC_content',
        'mappability': 'Mappability',
        'length': 'Length',
        'id': 'ID'
    })[OUTPUT_COLUMNS]

    if df_final.empty:
        raise ValueError(f'Dataset for sample {samplename} is empty. Check config file.')

    return df_final


def detect_genome_build(target_df: pd.DataFrame) -> bool:
    """
    Detect if using b37 genome build (no 'chr' prefix).

    Args:
        target_df: Target regions DataFrame

    Returns:
        True if b37 (no 'chr' prefix), False if GRCh38/hg38
    """
    # Check first data row (skip potential header)
    first_chrom = str(target_df.iloc[0, 0])
    if any(c.isalpha() for c in str(target_df.iloc[0, 2])):
        # First row is header, check second row
        first_chrom = str(target_df.iloc[1, 0])

    return not first_chrom.startswith('chr')


def main() -> int:
    """
    Main entry point for inference dataset generation.

    Returns:
        Exit code (0 for success, 1 for failure)
    """
    parser = argparse.ArgumentParser(
        description='Create annotated datasets for miXer inference'
    )
    parser.add_argument('-j', '--json', required=True,
                        help="Path to the miXer JSON config file")
    parser.add_argument('-s', '--samples', required=False,
                        help="Sample sheet file (optional, defaults to JSON config)")
    args = parser.parse_args()

    # Load config
    with open(args.json, 'r') as config_file:
        config = json.load(config_file)

    out_dir = os.path.join(
        os.path.abspath(config['main_outdir_host']),
        config['exp_id']
    )

    # Setup logging to output directory
    os.makedirs(out_dir, exist_ok=True)
    setup_logging(out_dir)

    # Validate config
    try:
        validate_config(config)
    except (ValueError, FileNotFoundError) as e:
        logging.error(f"Config validation failed: {e}")
        return 1

    # Determine samples file
    samples_file = args.samples if args.samples else config.get('sample_list')
    if not samples_file:
        logging.error("Samples file must be provided via -s argument or in JSON config")
        return 1

    # Load control NRC data
    pool_nrc_path = None
    if config.get('premade_control_rdata') and os.path.exists(config['premade_control_rdata']):
        pool_nrc_path = config['premade_control_rdata']
        logging.info(f"Using premade control file: {pool_nrc_path}")
    else:
        pool_nrc_path = os.path.join(
            out_dir, "excavator2_output", "output", "DataAnalysis_w50k",
            "Control", "RCNorm", "Control.NRC.RData"
        )
        logging.info(f"Using EXCAVATOR2 control file: {pool_nrc_path}")

    if not os.path.exists(pool_nrc_path):
        logging.error(f"Control RData file not found: {pool_nrc_path}")
        return 1

    pool_nrc = extract_nrc(pool_nrc_path)

    # Find sample RData files
    search_pattern = os.path.join(
        out_dir, "excavator2_output", "output", "DataPrepare_w50k",
        "*", "RCNorm", "*RData"
    )
    logging.info(f"Looking for RData files: {search_pattern}")
    sample_rdata_files = glob.glob(search_pattern)
    logging.info(f"Found {len(sample_rdata_files)} RData files")

    if not sample_rdata_files:
        logging.error(f"No RData files found in {search_pattern}")
        return 1

    # Load sample sheet and filter for test samples only
    samples_df = pd.read_table(samples_file, sep="\t")
    samples_df['sampleType'] = samples_df['sampleType'].str.strip().str.lower()
    samples_df['ID'] = samples_df['ID'].astype(str)

    # Only process test samples (sampleType == 't', exact match)
    test_samples = samples_df[samples_df['sampleType'] == 't']['ID'].tolist()
    logging.info(f"Found {len(test_samples)} test samples to process")

    if not test_samples:
        logging.warning("No test samples found in sample sheet. Nothing to process.")
        return 0

    # Load target file
    target_df = pd.read_csv(config['target'], sep="\t", index_col=False, header=None)
    target_df[0] = target_df[0].astype(str)

    # Check if using b37 (no 'chr' prefix)
    is_b37 = detect_genome_build(target_df)
    if is_b37:
        logging.info("Detected b37 genome build (no 'chr' prefix)")
    else:
        logging.info("Detected GRCh38/hg38 genome build ('chr' prefix)")

    # Create output directory
    date_str = datetime.now().strftime("%Y%m%d-%H%M%S")
    dataset_test_dir = os.path.join(out_dir, f'datasets_testing_{date_str}')
    os.makedirs(dataset_test_dir, exist_ok=True)

    # Create pybedtools temp directory
    tmp_folder = os.path.join(out_dir, f"{config['exp_id']}_miXer_tmp")
    os.makedirs(tmp_folder, exist_ok=True)
    helpers.set_tempdir(tmp_folder)

    # Prepare mappability data
    map_bt = contrib.bigwig.bigwig_to_bedgraph(config['map'])
    if is_b37:
        map_df = pd.read_table(map_bt.fn, names=['chr', 'start', 'end', 'mapp'])
        map_df['chr'] = map_df['chr'].str.replace('chr', '', regex=False)
        map_bt = BedTool.from_dataframe(map_df).sort()

    # Prepare gap and centromere data
    centro = pd.read_table(config['centro'])
    gap = pd.read_table(config['gap'])
    gap = gap[['chrom', 'chromStart', 'chromEnd']]
    centrogap = pd.concat([centro, gap]).sort_values(
        by=['chrom', 'chromStart', 'chromEnd']
    ).drop_duplicates().reset_index(drop=True)
    if is_b37:
        centrogap['chrom'] = centrogap['chrom'].str.replace('chr', '', regex=False)
    centrogap_bt = BedTool.from_dataframe(centrogap).sort()

    # Annotate target regions
    logging.info("Starting target annotation...")
    annotated_target = annotate_target(target_df, config['ref'], map_bt, centrogap_bt)

    # Pre-compute BedTool objects for reuse (performance optimization)
    pool_bed = BedTool.from_dataframe(pool_nrc).sort()
    target_bed = BedTool.from_dataframe(annotated_target).sort()

    # Process each test sample
    logging.info("Processing test samples...")
    processed_count = 0
    failed_samples = []

    for idx, rdata_file in enumerate(sample_rdata_files):
        sample_name = '.'.join(os.path.basename(rdata_file).split('.')[:-2])

        if sample_name not in test_samples:
            continue

        logging.info(f"Processing sample: {sample_name}")

        try:
            # Extract NRC and create dataset
            sample_nrc = extract_nrc(rdata_file)
            sample_dataset = make_dataset(
                sample_name, sample_nrc, pool_nrc, annotated_target,
                pool_bed=pool_bed, target_bed=target_bed
            )

            # Save dataset - single gzip-compressed write (fixes corruption bug)
            output_path = os.path.join(dataset_test_dir, f"{sample_name}_TARGET.txt.gz")
            sample_dataset.to_csv(output_path, sep='\t', index=False, compression='gzip')

            logging.info(f"Saved dataset: {output_path}")
            processed_count += 1

        except Exception as e:
            logging.error(f"Failed to process sample {sample_name}: {e}")
            failed_samples.append(sample_name)
            continue

        # Periodic cleanup every 10 samples to avoid temp file accumulation
        if (idx + 1) % 10 == 0:
            helpers.cleanup()
            gc.collect()

    # Final cleanup
    gc.collect()
    helpers.cleanup()
    shutil.rmtree(tmp_folder, ignore_errors=True)

    # Summary
    logging.info(f"Dataset generation complete. Processed: {processed_count}, Failed: {len(failed_samples)}")
    if failed_samples:
        logging.warning(f"Failed samples: {failed_samples}")

    return 0 if not failed_samples else 1


if __name__ == "__main__":
    sys.exit(main())
