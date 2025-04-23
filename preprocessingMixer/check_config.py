import pandas as pd
import argparse
import os
import sys

VALID_GENDERS = {"m", "f", "mf"}
#VALID_ANALYSIS = {"case", "ctrl"}
VALID_ANALYSIS = {"T", "C"}

parser = argparse.ArgumentParser(description="Validate and normalize config file for EXCAVATOR2 step")
parser.add_argument("-cf", "--config", required=True, help="Path to original config file")
parser.add_argument("-o", "--output_dir", required=True, help="Where to save the cleaned config")
parser.add_argument("-bd", "--bam_dir", required=False, help="Optional: base path to prepend to bam filenames")

args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)

# Load config
try:
    df = pd.read_csv(args.config, sep="\t", names=["ID", "bamPath", "Gender", "sampleType"], dtype=str)
except Exception as e:
    sys.exit(f"Failed to read config file: {e}")

# Normalize fields
df["ID"] = df["ID"].str.strip()
df["bamPath"] = df["bamPath"].str.strip()
df["Gender"] = df["Gender"].str.strip().str.lower()
df["analysis"] = df["analysis"].str.strip().str.lower()

# Validation: Gender and Analysis
invalid_gender = df[~df["Gender"].isin(VALID_GENDERS)]
invalid_analysis = df[~df["analysis"].isin(VALID_ANALYSIS)]

if not invalid_gender.empty:
    print("Invalid gender values found:")
    print(invalid_gender)

if not invalid_analysis.empty:
    print("Invalid analysis values found:")
    print(invalid_analysis)

if not invalid_gender.empty or not invalid_analysis.empty:
    sys.exit("Stopping due to invalid values in the config file.")

# Construct full BAM paths (if bam_dir provided)
if args.bam_dir:
    df["full_bam_path"] = df["bamPath"].apply(lambda x: os.path.join(args.bam_dir, x))
else:
    df["full_bam_path"] = df["bamPath"]

# Check BAM file existence
missing_bams = df[~df["full_bam_path"].apply(os.path.exists)]
if not missing_bams.empty:
    print("The following BAM files are missing:")
    print(missing_bams[["ID", "full_bam_path"]])
    sys.exit("Stopping due to missing BAM files.")

# Drop helper column before saving
df = df.drop(columns=["full_bam_path"])

# Save cleaned config
cleaned_name = os.path.basename(args.config) + "_forExca2"
output_path = os.path.join(args.output_dir, cleaned_name)
df.to_csv(output_path, sep="\t", index=False, header=True)

print(f"Cleaned config saved to: {output_path}")
