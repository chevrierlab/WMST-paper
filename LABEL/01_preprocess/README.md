

# Preprocesing Pipeline

This folder contains four scripts/notebooks for preparing data before model training and evaluation.

## 1. `01_label_file_preprocess.py`

- Merges all `*_annotation.tsv` files within each batch folder into a single CSV (`all_annotations.csv`).
- Default batches: `CTRL_1`, `CTRL_2`, `LPS_1`, `LPS_2`.

## 2. `02_TrainImageNormalization.ipynb`

- Learns image normalization parameters (e.g., stain normalization, mean–std statistics) from training images.
- Saves normalization settings for consistent preprocessing across datasets.

## 3. `03_DatasetGeneration.ipynb`

- Generates the main training/validation dataset.
- Applies normalization, crops/patches images, aligns labels, and creates dataset splits.

## 4. `04_AdditionalSetGeneration.ipynb`

- Creates additional evaluation datasets using the same preprocessing protocol.
- Ensures external or held-out samples follow the same normalization and formatting as the main dataset.


