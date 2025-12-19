#!/bin/bash
#SBATCH --job-name=ST
#SBATCH --time=30:50:00
#SBATCH --partition=caslake
#SBATCH --account=pi-nchevrier
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=100G
#SBATCH --output=output.out
#SBATCH --error=errorfile.err


# Path to STAR binary
export PATH=/project/nchevrier/project2_data/software/STAR_2.7.10a/STAR-2.7.10a/bin/Linux_x86_64_static/:$PATH

# Directory containing subfolders of FASTQ files
Project_Dir="../fastq"

Subfolders=("${Project_Dir}"/*/)
if [[ -d "${Subfolders[0]}" ]]; then
  echo "Detected multiple samples (subfolders). Proceeding with multi-sample alignment..."

    # Loop through each subfolder
    for Subfolder in "$Project_Dir"/*/; do
        # Get the name of the subfolder (basename removes the path and trailing slash)
        Subfolder_Name=$(basename "$Subfolder")

        # Define input FASTQ files for this subfolder
        Fastq_1=$(ls "${Subfolder}"*R1_*.fastq.gz 2>/dev/null | head -n 1)
        Fastq_2=$(ls "${Subfolder}"*R2_*.fastq.gz 2>/dev/null | head -n 1)

        # Check if FASTQ files exist
        if [[ -z "$Fastq_1" || -z "$Fastq_2" ]]; then
            echo "Missing files in $Subfolder. Skipping..."
            continue
        else
            echo "NO missing files in $Subfolder."
        fi
        
        # Create output directory named after the subfolder
        Output_Dir="ST/${Subfolder_Name}"
        mkdir -p "$Output_Dir"

        # Run STAR for this subfolder
        STAR \
        --runThreadN 24 \
        --genomeDir /project/nchevrier/projects/dcipurko/NGS_tools/refdata-gex-mm10-2020-A/star_2.7.10a \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist n12_barcodes.txt \
        --soloCBstart 1 \
        --soloCBlen 12 \
        --soloUMIstart 29 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloUMIdedup 1MM_CR \
        --soloCBmatchWLtype Exact \
        --soloFeatures Gene \
        --outFilterScoreMinOverLread 0.25 \
        --readFilesIn "$Fastq_2" "$Fastq_1" \
        --readFilesCommand zcat \
        --soloOutFileNames "${Output_Dir}/" features.tsv barcodes.tsv matrix.mtx

        # Compress the output
        gzip "${Output_Dir}/Gene/raw/"*

        # Clean up intermediate files if necessary
        rm Aligned.out.sam
    done
else
    # Define input FASTQ files directly
    Fastq_1=$(ls "${Project_Dir}"/*R1_*.fastq.gz 2>/dev/null | head -n 1)
    Fastq_2=$(ls "${Project_Dir}"/*R2_*.fastq.gz 2>/dev/null | head -n 1)

    # Check if FASTQ files exist
    if [[ -z "$Fastq_1" || -z "$Fastq_2" ]]; then
        echo "Missing FASTQ files in $Project_Dir. Exiting..."
        exit 1
    else
        echo "FASTQ files found: $Fastq_1 and $Fastq_2"
    fi

    # Define the output directory
    Output_Dir="ST"
    mkdir -p "$Output_Dir"

        # Run STAR for this subfolder
        STAR \
        --runThreadN 24 \
        --genomeDir /project/nchevrier/projects/dcipurko/NGS_tools/refdata-gex-mm10-2020-A/star_2.7.10a \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist n12_barcodes.txt \
        --soloCBstart 1 \
        --soloCBlen 12 \
        --soloUMIstart 29 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloUMIdedup 1MM_CR \
        --soloCBmatchWLtype Exact \
        --soloFeatures Gene \
        --outFilterScoreMinOverLread 0.25 \
        --readFilesIn "$Fastq_2" "$Fastq_1" \
        --readFilesCommand zcat \
        --soloOutFileNames "${Output_Dir}/" features.tsv barcodes.tsv matrix.mtx

        # Compress the output
        gzip "${Output_Dir}/Gene/raw/"*

        # Clean up intermediate files if necessary
        rm Aligned.out.sam
fi