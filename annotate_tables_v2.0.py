#!/usr/bin/env python3.11

"""
Python script that adds mapman annotations to differential expression result tables from nf-core/differentialabundance pipeline.

Will work for C.campestris and C.monogyna identifiers.
"""

try:
    import os
    import sys
    import pandas as pd
    import argparse
    from argparse import RawTextHelpFormatter
    from datetime import datetime
    import logging
except ModuleNotFoundError as e:
    print(f"Error: {e}")
    print("It seems that some required modules are not installed or the conda environment is not active.")
    print("Please ensure you have activated the appropriate conda environment and all dependencies are installed.")
    print("Use the following command to activate the correct environment:\n")
    print(f"    conda activate python_3.11\n")
    sys.exit(1)

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# Create logger
logger = logging.getLogger(__name__)



# Function to define species by checking string in the gene_id column. Will also set path for mapman annotation file based on species.
def check_species(input_file):
    with open(input_file, 'r') as file:
        header_line = file.readline()  # Read and discard the header line
        first_data_line = file.readline()  # Read the first line after the header
        if "Cmono" in first_data_line:
            species = "Cuscuta monogyna"
            mapman_path = '/Users/Shared/cuscuta/references/cuscuta_monogyna/cuscuta_monogyna.annot.cds.MapMan4_v7.0.txt'
        elif "Cc" in first_data_line:
            species = "Cuscuta campestris"
            mapman_path = '/Users/Shared/cuscuta/references/cuscuta_campestris/data_cucam_0.32.annot.cds.MapMan4_v7.0.txt'
        else:
            raise SystemExit("Species not recognized. Exiting!")
    return species, mapman_path


def process_mapman(mapman, col2):
    mapman = mapman.copy()

    # Capitalize first letter in gene identifiers (Mercator only outputs lower case for some reason)
    mapman[col2] = mapman[col2].str.replace(r'^.', lambda x: x.group().upper(), regex = True)
    # Remove rows with no identifier
    mapman = mapman.dropna(subset=[col2])
    # Extract iso number
    mapman['iso'] = mapman[col2].str.extract(r'(?:\.t|\.)(\d+)$').astype(int)
    # Remove .t1, .t2, etc. from IDENTIFIER
    mapman[col2] = mapman[col2].str.replace(r'(?:\.t|\.)(\d+)$', '', regex=True)
    # Group by IDENTIFIER and sort by iso number
    sorted_mapman = mapman.sort_values(by='iso').groupby(col2)

    # Create a new DataFrame to store the results
    result = pd.DataFrame(columns=mapman.columns)
    # Iterate over unique identifiers and keep rows with the lowest iso number
    for identifier, group in sorted_mapman:
        min_iso = group['iso'].min()
        result = pd.concat([result, group[group['iso'] == min_iso]])

    result = result.reset_index(drop=True)
    return result

def annotate_table(dataset, col1, mapman, col2):
    # Check if all identifiers in dataset are also in the annotation table
    missing_annotations = dataset[~dataset[col1].isin(mapman[col2])]
    if not missing_annotations.empty:
        logger.info("Some identifiers are missing annotation.")

    # Keep annotations for identifiers in the dataset
    annotation = mapman[mapman[col2].isin(dataset[col1])]
    # Record number of bins in annotation per identifier
    bins = annotation[col2].value_counts()
    bins = bins.sort_index()
    
    # Record number of occurrences in dataset per identifier
    freq = dataset[col1].value_counts()
    freq = freq.sort_index()

    # Create replication vector for dataset
    rep_dataset = bins.reindex(freq.index.repeat(freq)).fillna(0).astype(int)
    dataset = dataset.loc[dataset.index.repeat(rep_dataset)]
    
    # Create replication vector for annotation
    rep_annotation = freq.reindex(bins.index.repeat(bins)).fillna(0).astype(int)
    annotation = annotation.loc[annotation.index.repeat(rep_annotation)]
    
    # Ensure same order
    dataset = dataset.sort_values(by=col1)
    annotation = annotation.sort_values(by=col2)

    # Merge dataset and annotation table
    annotated_dataset = pd.concat([dataset.reset_index(drop=True), annotation.reset_index(drop=True)], axis=1)
    # Redefine factor levels
    annotated_dataset[col1] = pd.Categorical(annotated_dataset[col1], categories=sorted(annotated_dataset[col1].unique()))
    annotated_dataset[col2] = pd.Categorical(annotated_dataset[col2], categories=sorted(annotated_dataset[col2].unique()))

    # Test if all identifiers on a row are the same
    if not annotated_dataset[col1].equals(annotated_dataset[col2]):
        logger.info(f"Mismatches detected!: Exiting!")
        sys.exit(1)
    else:
        logger.info(f"No mismatches found!")

    
        

    # Return annotated dataset
    return annotated_dataset


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Annotate a table based on a reference table.', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-indir', type = str, default = '03_differentialabundance_out/tables/differential', help = 'Path to directory with DEG tables')
    parser.add_argument('-col1', type=str, default = 'gene_id', help = '''Column name with identifiers in the dataset.
    Default: gene_id''')
    #parser.add_argument('-mapman', type=str, help = 'Path to MapMan annotation reference table')
    parser.add_argument('-col2', type=str, default = 'IDENTIFIER', help = '''Column name with identifiers in the reference table.
    Default: IDENTIFIER''')
    parser.add_argument('-outdir', type=str, default = '04_annotated_DE_tables', help = 'Output directory')
    parser.add_argument('-lfc', type=float, default = 1.5, help = '''LogFoldChange cutoff for differentially expressed genes.
    Default: 1.5''')
    parser.add_argument('-padj', type=float, default = 0.05, help = '''Adjusted p-value cutoff,
    Default: 0.05''')
    args = parser.parse_args()

    # Log startup info
    logger.info("Starting the MapMan annotation script.")
    logger.info(f"Input directory: {args.indir}")
    logger.info(f"Output directory: {args.outdir}")
    logger.info(f"Column with dataset identifiers: {args.col1}")
    logger.info(f"Column with reference identifiers: {args.col2}")
    logger.info(f"LogFoldChange cutoff: {args.lfc}")
    logger.info(f"Adjusted p-value cutoff: {args.padj}")

    # Check gene ids of the first file in the input directory to check species
    files = os.listdir(args.indir)    
    # Ensure there is at least one file in the directory
    if not files:
        raise SystemExit("No files found in the specified directory.")
    # Select the first file to check the species
    first_file = files[0]
    # Check species (read in first file in input directory and check species)
    species, mapman_path = check_species(os.path.join(args.indir, first_file))

    logger.info(f"Checking input table")
    logger.info(f"Using annotation file for {species}")
   

    # Check a line from mapman table to determine quote character of the file
    with open(mapman_path, 'r') as file:
        first_line = file.readline()
        second_line = file.readline()
        if '"' in second_line:
            quote_char = '"'
        elif "'" in second_line:
            quote_char = "'"
        else:
            quote_char = '"'


    # Processed mapman annotation file
    mapman = pd.read_csv(mapman_path, delimiter='\t', quotechar = quote_char)
    mapman_filename = os.path.basename(mapman_path)
    logger.info(f"Processing MapMan annotation file: {mapman_filename}")
    mapman = process_mapman(mapman, args.col2)

    # Read in input differential expression files and add mapman annotations
    input_dir = args.indir
    for filename in os.listdir(input_dir):
        file_path = os.path.join(input_dir, filename)

        if os.path.isfile(file_path):
            output_basename = os.path.splitext(filename)[0]
            logger.info(f"Adding MapMan annotations to dataset:\t{filename}")
            dataset = pd.read_csv(file_path, delimiter = '\t')
            annotated_data = annotate_table(dataset, args.col1, mapman, args.col2)
            # Remove the two last columns
            annotated_data = annotated_data.drop(columns = ["TYPE", "iso"])

            # Filter significant genes (log2FC > 2 & padj < 0.05)
            annotated_data_sig = annotated_data[(abs(annotated_data["log2FoldChange"]) > args.lfc) & (annotated_data["padj"] < args.padj)]

            # Save annotated dataset to a tab delimited file and excel file
            logger.info(f"Saving output tables in {args.outdir}")
            logger.info(f"Writing tsv: {output_basename}.annot.tsv")
            logger.info(f"Writing xlsx: {output_basename}.annot.xlsx")
            annotated_data.to_csv(args.outdir + '/' + output_basename + '.annot.tsv', sep='\t', index=False)
            annotated_data.to_excel(args.outdir + '/' + output_basename + '.annot.xlsx', index = False)
    
            # Save subset of significant genes to tab delimited and excel file
            annotated_data_sig.to_csv(args.outdir + '/' + output_basename + '.annot.sig.tsv', sep='\t', index=False)
            annotated_data_sig.to_excel(args.outdir + '/' + output_basename + '.annot.sig.xlsx', index = False)
            logger.info(f"Subsetting significantly differentially expressed genes")
            logger.info(f"Number of significantly differentially expressed genes: {len(annotated_data_sig)}")
            logger.info(f"Writing tsv: {output_basename}.annot.sig.tsv")
            logger.info(f"Writing xlsx: {output_basename}.annot.sig.xlsx")
    
    logger.info("FINISHED!")
