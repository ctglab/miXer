import os
import argparse
import json
import glob
import sys
import pyreadr
import numpy as np
import pandas as pd
import logging

import pyBigWig
from pybedtools import BedTool

from typing import Dict, Tuple

# setup logger
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler('mixerCreateDatasets.log'),
        logging.StreamHandler(sys.stdout)
    ])

# Log uncaught exceptions
def handle_exception(exc_type, exc_value, exc_traceback):
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))

sys.excepthook = handle_exception

#### fix_target: keep only chromosomal coordinates (remove additional columns if present)
def fix_target(a, min_TR_size=10):
    df = a.copy()
    
    # Identify and handle header row
    start_row = 0
    if any(c.isalpha() for c in str(df.iloc[0, 2])):
        start_row = 1
    
    # Coordinates are expected in columns 1 and 2
    data = df.iloc[start_row:].copy()
    data_coords = data.iloc[:, 1:3].astype(int)
    mask = (data_coords.iloc[:, 1] - data_coords.iloc[:, 0]) >= min_TR_size
    
    # Keep only the filtered data rows
    target_slice = data[mask].reset_index(drop=True)
    
    # Rename coordinates columns for consistency with downstream logic (BedTool needs names or specific structure)
    new_cols = ['chrom', 'start', 'end'] + list(target_slice.columns[3:])
    target_slice.columns = new_cols

    target_bed = BedTool.from_dataframe(target_slice).sort()
    if not target_bed or target_bed == "0":
       raise ValueError(f"Check that target file is properly formatted and contains regions >= {min_TR_size} bp")
    else:
        return target_bed

def to_bw_chrom(ch: str, target_has_chr: bool, bw_has_chr: bool) -> str:
    """
    Normalize chromosome naming between the target (df rows) and bigWig.
    - If target uses 'chr' and bigWig doesn't → drop 'chr' before querying BW.
    - If target doesn't use 'chr' and bigWig does → add 'chr' before querying BW.
    - Else return as-is.
    """
    ch = str(ch)
    if target_has_chr and not bw_has_chr:
        return ch[3:] if ch.startswith("chr") else ch
    if not target_has_chr and bw_has_chr:
        return ch if ch.startswith("chr") else "chr" + ch
    return ch

def normalize_chr_series(s: pd.Series, want_chr: bool) -> pd.Series:
    """Normalize chromosome names to consistently use or omit the 'chr' prefix.
    Handles mixed inputs (e.g. ['chr1', '2', 'chrX']) correctly per-element."""
    s = s.astype(str)
    stripped = s.str.replace(r'^chr', '', n=1, regex=True)
    return 'chr' + stripped if want_chr else stripped


def _read_bedtool_result(bedtool_result, columns, dtypes=None):
    """Safely read a BedTool result file, returning an empty DataFrame
    with the expected columns if the intersection produced no results."""
    if dtypes is None:
        dtypes = {}
    try:
        df = pd.read_table(bedtool_result.fn, header=None, dtype=dtypes)
        if df.empty:
            return pd.DataFrame(columns=columns)
        df.columns = columns
        return df
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)

def annotate_target(tar, ref, mapp_bw_path, gapcen, min_TR_size=10):
       
    ## GC content
    nuc = BedTool.nucleotide_content(fix_target(tar, min_TR_size), fi=ref)
    df_with_gc = pd.read_table(nuc.fn, dtype={0: str, 1: int, 2: int}).loc[:, ['#1_usercol', '2_usercol', '3_usercol', '5_pct_gc']]
    df_with_gc.columns = ['chrom', 'start', 'end', 'gc']

    # Mappability using pyBigWig
    bw = pyBigWig.open(mapp_bw_path)
    bw_chroms = set(bw.chroms().keys())

    # Detect naming conventions
    target_has_chr = df_with_gc['chrom'].astype(str).str.startswith('chr').any()
    bw_has_chr = any(c.startswith('chr') for c in bw_chroms)

    mappabilities = []
    for _, row in df_with_gc.iterrows():
        bw_ch = to_bw_chrom(row['chrom'], target_has_chr, bw_has_chr)
        try:
            if bw_ch in bw_chroms:
                vals = bw.values(bw_ch, int(row['start']), int(row['end']))
                # remove None and NaN
                vals = [v for v in vals if (v is not None and not np.isnan(v))]
                mean_val = float(sum(vals) / len(vals)) if vals else 0.0
            else:
                mean_val = 0.0
        except Exception:
            mean_val = 0.0
        mappabilities.append(mean_val)

    bw.close()

    df_with_gc["mappability"] = mappabilities

    ## Make BedTool to remove GAP regions
    df_with_mapp_bed = BedTool.from_dataframe(df_with_gc).sort()
    final_target = df_with_mapp_bed.intersect(gapcen, v=True)
    final_target_df = pd.read_table(
        final_target.fn,
        names=['Chr', 'Start', 'End', 'GC_content', 'Mappability'], dtype={0: str, 1: int, 2: int}
    )

    final_target_df['Length'] = final_target_df['End'] - final_target_df['Start']

    if final_target_df.empty:
        raise ValueError("Check that target file is not empty and properly formatted")
    return final_target_df

def make_dataset(sample_id: str,
                 sample_df: pd.DataFrame,
                 control_pool_df: pd.DataFrame,
                 target_annot_df: pd.DataFrame,
                 drop_offtarget: bool = True) -> pd.DataFrame:
    """
    Create a normalized dataset by aligning sample NRC values with target regions
    and normalizing using control pool NRCs. Keeps both on-target and off-target regions.

    Parameters:
        sample_id (str): Sample identifier.
        sample_df (pd.DataFrame): Must contain ['chrom', 'start', 'end', 'RCNorm', 'INOUT', 'source_file'].
        control_pool_df (pd.DataFrame): Must contain ['chrom', 'start', 'end', 'RCNorm_pool', 'INOUT'].
        target_annot_df (pd.DataFrame): Must contain ['Chr', 'Start', 'End', 'GC_content', 'Mappability', 'Length'].

    Returns:
        pd.DataFrame: Annotated, normalized DataFrame with region classification.
    """

    # Ensure correct dtypes
    for df in [sample_df, control_pool_df, target_annot_df]:
        for col in ['start', 'end']:
            if col in df.columns:
                df[col] = df[col].astype(int)
        if 'chrom' in df.columns:
            df['chrom'] = df['chrom'].astype(str)

    # Step 1: Intersect sample and pool with targets
    sample_bt = BedTool.from_dataframe(sample_df[['chrom', 'start', 'end', 'RCNorm']]).sort()
    pool_bt = BedTool.from_dataframe(control_pool_df[['chrom', 'start', 'end', 'RCNorm_pool']]).sort()
    target_bt = BedTool.from_dataframe(target_annot_df[['Chr', 'Start', 'End']]).sort()

    # Sample on/off-target
    intersect_sample_on = sample_bt.intersect(target_bt, wo=True)
    intersect_sample_off = sample_bt.intersect(target_bt, v=True)

    sample_on = pd.read_table(intersect_sample_on.fn, header=None, dtype={0: str, 1: int, 2: int}) # ensure chrom is string, coords are int
    sample_on.columns = ['chrom_sample', 'start_sample', 'end_sample', 'RCNorm', 'chrom_tar', 'start_tar', 'end_tar', 'overlap']
    sample_on = sample_on.rename(columns={'chrom_tar': 'chrom', 'start_tar': 'start', 'end_tar': 'end'})
    sample_on = sample_on[['chrom', 'start', 'end', 'RCNorm']]

    sample_off = pd.read_table(intersect_sample_off.fn, header=None, dtype={0: str, 1: int, 2: int})
    sample_off.columns = ['chrom', 'start', 'end', 'RCNorm']

    sample_annot = pd.concat([sample_on, sample_off], ignore_index=True)

    # Pool on/off-target
    intersect_pool_on = pool_bt.intersect(target_bt, wo=True)
    intersect_pool_off = pool_bt.intersect(target_bt, v=True)

    # ensure chrom is string, coords are int
    pool_on = pd.read_table(intersect_pool_on.fn, header=None, dtype={0: str, 1: int, 2: int})
    pool_on.columns = ['chrom_sample', 'start_sample', 'end_sample', 'RCNorm_pool', 'chrom_tar', 'start_tar', 'end_tar', 'overlap']
    pool_on = pool_on.rename(columns={'chrom_tar': 'chrom', 'start_tar': 'start', 'end_tar': 'end'})
    pool_on = pool_on[['chrom', 'start', 'end', 'RCNorm_pool']]

    pool_off = pd.read_table(intersect_pool_off.fn, header=None, dtype={0: str, 1: int, 2: int})
    pool_off.columns = ['chrom', 'start', 'end', 'RCNorm_pool']
    pool_annot = pd.concat([pool_on, pool_off], ignore_index=True)

    # Step 2: Merge sample and control
    merged_df = pd.merge(sample_annot, pool_annot, on=['chrom', 'start', 'end'])

    # Step 3: Merge with target annotation (on-target only will match)
    merged_df = pd.merge(
        merged_df,
        target_annot_df.rename(columns={'Chr': 'chrom', 'Start': 'start', 'End': 'end'}),
        on=['chrom', 'start', 'end'],
        how='left'  # allows off-target rows to remain
    )
    # Step 4: Normalize and log-transform
    merged_df['NRC_poolNorm'] = merged_df['RCNorm'] / merged_df['RCNorm_pool']
    merged_df['NRC_poolNorm'] = np.log2(merged_df['NRC_poolNorm'].replace(0, np.nan))

    merged_df['id'] = sample_id

    # Step 5: Add sample metadata
    sample_df['start'] = sample_df['start'].astype(int)
    sample_df['end'] = sample_df['end'].astype(int)
    sample_df['chrom'] = sample_df['chrom'].astype(str)

    sample_info = sample_df[['chrom', 'start', 'end', 'INOUT', 'source_file']]
    final_df = pd.merge(merged_df, sample_info, on=['chrom', 'start', 'end'], how='left')

    # Step 6: Final columns
    final_df = final_df[[
        'chrom', 'start', 'end', 'GC_content', 'Mappability', 'Length',
        'NRC_poolNorm', 'INOUT'
    ]]

    # Normalize the median (don't consider off target and X/Y chrs)
    nrc_median = final_df[(~final_df["chrom"].isin(["chrX", "ChrX", "chrx", "X", 
                                                    "chrY", "ChrY", "chry", "Y"])) &
                                                    (final_df["INOUT"] != "OUT")]["NRC_poolNorm"].median()
    
    # Normalize the median (don't consider off target and X/Y chrs)
    nrc_median_off = final_df[(~final_df["chrom"].isin(["chrX", "ChrX", "chrx", "X", 
                                                    "chrY", "ChrY", "chry", "Y"])) &
                                                    (final_df["INOUT"] == "OUT")]["NRC_poolNorm"].median()
    

    #final_df['NRC_poolNorm'] = final_df['NRC_poolNorm'] - nrc_median
    final_df["NRC_poolNorm_autosomalMedian"] = nrc_median
    final_df["NRC_poolNorm_offtarget_autosomalMedian"] = nrc_median_off
    # drop offtarget from final_df
    if drop_offtarget:
        final_df = final_df.loc[final_df['INOUT'] == 'IN']
    return final_df, nrc_median, nrc_median_off

def normalize_by_autosomal_median(all_dataframes: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, pd.DataFrame], Dict[str, str]]:
    """
    Normalizes samples using the autosomal median for on-target (INOUT=="IN")
    and off-target (INOUT=="OUT") positions.

    Parameters:
        all_dataframes (dict): map sample_id -> DataFrame. Each DataFrame
            must have at least these columns:
              - 'NRC_poolNorm'                     (float)
              - 'INOUT'                            (string: "IN" or "OUT")
              - 'NRC_poolNorm_autosomalMedian'     (float, the on-target median)
              - 'NRC_poolNorm_offtarget_autosomalMedian' (float, the off-target median)

    Returns:
        normalized_dataframes (dict): map sample_id -> DataFrame with additional
            column "NRC_poolNorm_normalized" (float).
        quarantined_samples (dict): map sample_id -> quarantine reason
            (currently unused in this simplified version).
    """
    normalized_dataframes: Dict[str, pd.DataFrame] = {}
    quarantined_samples: Dict[str, str]   = {}

    # 1) Gather sample medians
    medians_dict: Dict[str, Dict[str, float]] = {}
    for sample_id, df in all_dataframes.items():
        if "NRC_poolNorm_autosomalMedian" not in df.columns or \
           "NRC_poolNorm_offtarget_autosomalMedian" not in df.columns:
            # If the required columns are missing, we cannot normalize
            continue

        # Get the medians from the first row (they are constant for the entire DataFrame)
        med_on  = df["NRC_poolNorm_autosomalMedian"].iloc[0]
        med_off = df["NRC_poolNorm_offtarget_autosomalMedian"].iloc[0]
        medians_dict[sample_id] = {"on": med_on, "off": med_off}

    logging.info("Starting samples normalization")

    # 2) For each sample, apply normalization
    for sample_id, df in all_dataframes.items():
        # If the required medians are missing, skip
        if sample_id not in medians_dict:
            logging.warning(f"Median values for sample {sample_id} not found, skipping.")
            continue

        df_copy = df.copy()

        # Get the medians for this sample
        med_on  = medians_dict[sample_id]["on"]
        med_off = medians_dict[sample_id]["off"]

        logging.info(f"Applying autosomal median normalization to sample {sample_id}: Current Autosomal median IN-target = {med_on:.4f}, OUT-target = {med_off:.4f}")

        # Create masks
        mask_in  = df_copy["INOUT"] == "IN"
        mask_out = df_copy["INOUT"] == "OUT"

        df_copy.loc[mask_in,  "NRC_poolNorm_normalized"] = df_copy.loc[mask_in,  "NRC_poolNorm"] - med_on
        df_copy.loc[mask_out, "NRC_poolNorm_normalized"] = df_copy.loc[mask_out, "NRC_poolNorm"] - med_off
        
        # Print additional statistics
        auto_in_median = df_copy[(df_copy["INOUT"] == "IN") & (~df_copy["chrom"].str.upper().isin(["CHRX", "X"]))]["NRC_poolNorm_normalized"].median()
        auto_out_median = df_copy[(df_copy["INOUT"] == "OUT") & (~df_copy["chrom"].str.upper().isin(["CHRX", "X"]))]["NRC_poolNorm_normalized"].median()
        chrX_median = df_copy[df_copy["chrom"].str.upper().isin(["CHRX", "X"])]["NRC_poolNorm_normalized"].median()

        logging.info(f"Post-normalization medians (Autosomes IN TRs | Autosomes OUT TRs | chrX TRs): {auto_in_median:.4f} / {auto_out_median:.4f} / {chrX_median:.4f}")

        normalized_dataframes[sample_id] = df_copy

    return normalized_dataframes, quarantined_samples

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create annotated datasets for miXer')
    parser.add_argument('-j', '--json', help="Path to the miXer json file", required=True)

    args = parser.parse_args()
    
    # Load configuration
    with open(args.json, 'r') as j:
        config = json.load(j)
    
    # 1. Setup Directories
    outdir_host = config.get('main_outdir_host')
    exp_id = config.get('exp_id')
    
    if not outdir_host or not exp_id:
        raise ValueError("Missing 'main_outdir_host' or 'exp_id' in JSON config.")

    # Base output directory
    base_outdir = os.path.join(os.path.abspath(outdir_host), exp_id)
    if not os.path.exists(base_outdir):
        os.makedirs(base_outdir)
        
    dataset_outdir = os.path.join(base_outdir, f"datasets_miXer")
    os.makedirs(dataset_outdir, exist_ok=True)
    
    logging.info(f"Output directory: {dataset_outdir}")

    # 2. Load Support Files from JSON
    # Reference
    genome_fasta = config.get('ref')
    if not genome_fasta or not os.path.exists(genome_fasta):
        raise FileNotFoundError(f"Reference FASTA not found: {genome_fasta}")
        
    # Target
    target_file = config.get('target')
    if not target_file or not os.path.exists(target_file):
        raise FileNotFoundError(f"Target file not found: {target_file}")
    
    # Check for filtered target version
    min_TR_size = config.get('min_TR_size', 10)
    filtered_target_dir = os.path.join(base_outdir, f"filtered_target_TR_geq_{min_TR_size}")
    filtered_target_path = os.path.join(filtered_target_dir, os.path.basename(target_file))
    
    if os.path.exists(filtered_target_path):
        logging.info(f"Using filtered target file found at: {filtered_target_path}")
        target_file = filtered_target_path
    else:
        logging.info(f"Filtered target file not found at {filtered_target_path}. Using original target: {target_file}")
        
    # Mappability
    map_bw = config.get('map')
    if not map_bw or not os.path.exists(map_bw):
        raise FileNotFoundError(f"Mappability BigWig not found: {map_bw}")

    # Centromeres & Gaps
    centro_file = config.get('centro')
    gap_file = config.get('gap')
    if not centro_file or not os.path.exists(centro_file):
        raise FileNotFoundError(f"Centromere file not found: {centro_file}")
    if not gap_file or not os.path.exists(gap_file):
         raise FileNotFoundError(f"Gap file not found: {gap_file}")

    # 3. Load Sample List
    sample_list_path = config.get('sample_list')
    if not sample_list_path or not os.path.exists(sample_list_path):
        raise FileNotFoundError(f"Sample list not found: {sample_list_path}")

    # Reading sample list assuming standard format: IDs, bamPath, Gender, sampleType
    scan_df = pd.read_csv(sample_list_path, sep="\t", dtype=str)

    logging.info(f"Loaded {len(scan_df)} samples from {sample_list_path}")

    # 4. Prepare Support Data
    logging.info("Loading support data...")
    # Target
    target = pd.read_csv(target_file, sep="\t", header=None, comment="#", usecols=[0, 1, 2], names=["chr", "start", "end"])
    target_has_chr = target['chr'].astype(str).str.startswith('chr').any()
    
    # Centromere & Gap
    centro = pd.read_csv(centro_file, sep="\t") 
    # check columns for centro -> expects chrom, chromStart, chromEnd
    if "chrom" not in centro.columns: # fallback if standard bed-like or different header
         centro = pd.read_csv(centro_file, sep="\t", header=None, names=["chrom", "chromStart", "chromEnd"])

    gap = pd.read_csv(gap_file, sep="\t")
    if "chrom" not in gap.columns:
        gap = pd.read_csv(gap_file, sep="\t", header=None, names=["chrom", "chromStart", "chromEnd"])
    gap = gap[['chrom', 'chromStart', 'chromEnd']]
    
    centrogap = pd.concat([centro, gap]).sort_values(by=['chrom', 'chromStart', 'chromEnd']).drop_duplicates().reset_index(drop=True)
    checkb37 = target[target['chr'].str.startswith(('chr'))]
    if checkb37.empty:
        centrogap['chrom'] = centrogap['chrom'].str.replace('chr','')
    
    centrogap_bt = BedTool.from_dataframe(centrogap).sort()

    # 5. Annotate Target
    min_TR_size = config.get('min_TR_size', 10)
    logging.info(f"Annotating target (min_TR_size={min_TR_size})...")
    annotated_target = annotate_target(target, genome_fasta, map_bw, centrogap_bt, min_TR_size)
    logging.info("Target annotation complete.")

    # 6. Locate Input RData Files
    # Path logic: <outdir>/<exp_id>/excavator2_output/output/DataPrepare_w50k/*/RCNorm/*RData
    rdata_search_pattern = os.path.join(base_outdir, "excavator2_output", "output", "DataPrepare_w50k", "*", "RCNorm", "*.RData")
    logging.info(f"Searching for RData files in: {rdata_search_pattern}")
    rdata_files = glob.glob(rdata_search_pattern)
    
    if not rdata_files:
        logging.error("No RData files found. Check preprocessing output.")
        sys.exit(1)
        
    t_count = 0
    c_count = 0
    for filepath in rdata_files:
        fname_base = os.path.basename(filepath)
        for _, row in scan_df.iterrows():
            sim_id = row['ID']
            if sim_id in fname_base:
                if fname_base == sim_id + ".RData" or fname_base.startswith(sim_id + "."):
                    stype = str(row['sampleType']).lower()
                    if stype == 't':
                        t_count += 1
                    elif stype == 'c':
                        c_count += 1
                    break
    
    logging.info(f"Found {len(rdata_files)} RData files ({t_count} tagged as 't', {c_count} tagged as 'c').")
    
    # 7. Locate Control Pool
    # Logic: explicit in config OR default path
    pool_rdata_path = config.get('premade_control_rdata')
    if pool_rdata_path:
        logging.info(f"Using external premade control RData provided in config: {pool_rdata_path}")
    else:
        # Default path
        logging.info("No premade control RData provided. Using default path.")
        pool_rdata_path = os.path.join(base_outdir, "excavator2_output", "output", "DataAnalysis_w50k", "Control", "RCNorm", "Control.NRC.RData")
    
    if not os.path.exists(pool_rdata_path):
        logging.error(f"Control RData not found at: {pool_rdata_path}")
        sys.exit(1)
        
    # Read Control Data
    try:
        controlPool_data_raw = pyreadr.read_r(pool_rdata_path)
        controlPool_data = next(iter(controlPool_data_raw.values()))
        controlPool_data = controlPool_data.rename(columns={'Class': 'INOUT', 'RCNorm': 'RCNorm_pool'})
        controlPool_data['chrom'] = normalize_chr_series(controlPool_data['chrom'], want_chr=target_has_chr)
    except Exception as e:
        raise RuntimeError(f"Error reading control RData: {e}") from e

    # 8. Process Samples
    
    # Iterate over found RData files
    for filepath in rdata_files:
        filename = os.path.basename(filepath)
        # Assuming filename is like SampleID.something.RData or just SampleID.RData
        # Does filename start with ID?
        
        file_id = None
        # Try exact match strategies (assuming standard excavator naming conventions)
        # Strategy A: filename stem (remove extension)
        # Strategy B: more aggressive separate by dot
        
        # Let's look up in the scan_df
        found_id_in_df = False
        
        # Helper to find ID
        fname_base = os.path.basename(filepath)
        
        for sim_id in scan_df['ID']:
            # Exact match in filename parts?
            # if sim_id is "Sample1" and filename is "Sample1.RData" or "Sample1.NRC.RData"
            if sim_id in fname_base:
                # Check if it's a "clean" match (e.g. followed by dot or end)
                # This avoids matching "Sample1" in "Sample10.RData"
                if fname_base == sim_id + ".RData" or fname_base.startswith(sim_id + "."):
                    file_id = sim_id
                    found_id_in_df = True
                    break
        
        if not found_id_in_df:
            logging.warning(f"Skipping {fname_base}: ID not found in sample list.")
            continue
            
        # Get metadata
        row = scan_df[scan_df['ID'] == file_id].iloc[0]
        gender = row['Gender']
        stype = str(row['sampleType']).lower()
        
        if stype == 'c': #skip control samples
            continue

        logging.info(f"Processing sample {file_id}")
        
        # Read Result
        try:
            res = pyreadr.read_r(filepath)
            df_samp = next(iter(res.values()))
            df_samp['chrom'] = normalize_chr_series(df_samp['chrom'], want_chr=target_has_chr)
            df_samp['source_file'] = file_id
            df_samp = df_samp.rename(columns={'Class': 'INOUT'})
        except Exception as e:
            logging.error(f"Error reading {filepath}: {e}")
            continue

        # Make Dataset
        try:
            mixer_data, nrc_aut_med, nrc_off_med = make_dataset(
                sample_id=file_id, 
                sample_df=df_samp, 
                control_pool_df=controlPool_data, 
                target_annot_df=annotated_target
            )
        except Exception as e:
            logging.error(f"Error creating dataset for {file_id}: {e}")
            continue

        # Local Normalization for single test sample
        one_div = {file_id: mixer_data.copy()}
        norm_dict, quar = normalize_by_autosomal_median(one_div)
        
        final_df = norm_dict.get(file_id, mixer_data.copy())
        if file_id not in norm_dict:
            final_df["NRC_poolNorm_normalized"] = np.nan
        
        # rename columns to match R notation
        final_df.rename(
            columns={'chrom': 'Chr',
            'start': 'Start',
            'end': 'End'},
            inplace=True)
        final_df.to_csv(os.path.join(dataset_outdir, f"{file_id}_miXer_data.tsv"), sep="\t", index=False)
        

    logging.info("Processing complete.")

