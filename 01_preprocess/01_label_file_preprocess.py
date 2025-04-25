import pandas as pd
import glob

# Define the list of batches to process
batch_list = ['CTRL_1', 'CTRL_2', 'LPS_1', 'LPS_2']

for batch in batch_list:
    # Find all TSV annotation files within the batch folder
    tsv_files = glob.glob(f'../Celltype_Annotations/{batch}/*_annotation.tsv')

    # Initialize an empty DataFrame to hold merged data
    all_data = pd.DataFrame()

    # Read and concatenate each TSV file
    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep='\t')  # Read using tab separator
        all_data = pd.concat([all_data, df], ignore_index=True)

    # Save the combined annotations to a single CSV file for the batch
    all_data.to_csv(f'../Celltype_Annotations/{batch}/all_annotations.csv', index=False)

    print("Done")

