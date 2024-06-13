import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Convert .bed file to .parquet format")
    parser.add_argument("--input_bed", required=True, help="Path to the input .bed file")
    parser.add_argument("--output_parquet", required=True, help="Path to the output .parquet file")
    return parser.parse_args()

def main():
    args = parse_args()
    
    # basic stats
    column_names = ["chrom", "start", "end", "snp", "score", "relation"]

    df = pd.read_csv(args.input_bed, sep='\t', header=None, names=column_names)

    df.to_parquet(args.output_parquet, index=False)

if __name__ == "__main__":
    main()
