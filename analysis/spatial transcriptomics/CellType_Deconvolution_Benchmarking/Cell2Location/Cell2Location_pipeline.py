import os
import numpy as np
import pandas as pd
import anndata as ad
ad.settings.allow_write_nullable_strings = True
from scipy import io
import sys
import types
import torch, scvi
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel, Cell2location
import cell2location
from scipy import sparse as sp
import os
import itertools
import pandas as pd
import json
import numpy as np
import pandas as pd



# ----------------- Config -----------------
os.environ["OMP_NUM_THREADS"] = "10"
os.environ["OPENBLAS_NUM_THREADS"] = "10"
os.environ["MKL_NUM_THREADS"] = "10"
torch.set_num_threads(10)
torch.set_num_interop_threads(2)
scvi.settings.dl_num_workers = 2 
scvi.settings.seed = 0

organs = ["SP", "LU", "LI", "KI"] 
samples = ["CTRL_1", "CTRL_2", "LPS_1", "LPS_2"]
combos = list(itertools.product(organs, samples))


array_id = int(os.environ.get("ARRAY_ID", "0"))
organ, sample = combos[array_id]
print(f"Selected Organ: {organ}, Sample: {sample}")

base_dir = "/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Reference_data"
comp_dir = os.path.join(base_dir, organ, "components")
out_dir  = os.path.join("/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/cell2location_output", organ)
os.makedirs(out_dir, exist_ok=True)




SAFE_SCALARS = (str, int, float, bool, np.generic)

def _stringify_complex(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    if isinstance(x, SAFE_SCALARS):
        return str(x)
    # Lists, dicts, sets, arrays, tuples, custom objects -> JSON string
    try:
        if isinstance(x, (list, dict, set, tuple, np.ndarray)):
            return json.dumps(x, default=str, ensure_ascii=False)
        # Last resort
        return json.dumps(x, default=str, ensure_ascii=False)
    except Exception:
        return str(x)

def sanitize_obs_for_h5ad(adata, verbose=True):
    """
    Coerce all obs columns to h5py-writable dtypes.
    - Categorical -> string
    - Datetime/Timedelta -> ISO string
    - Object columns are stringified elementwise if they contain any non-scalar values
    - Ensure obs_names are strings
    """
    obs = adata.obs.copy()

    # Ensure index are plain strings
    adata.obs_names = adata.obs_names.astype(str)

    converted_cols = []
    for col in obs.columns:
        s = obs[col]

        # 1) Categorical
        if pd.api.types.is_categorical_dtype(s):
            obs[col] = s.astype(str)
            converted_cols.append((col, "category→string"))
            continue

        # 2) Datetime-like
        if pd.api.types.is_datetime64_any_dtype(s) or pd.api.types.is_timedelta64_dtype(s):
            obs[col] = s.astype("string").fillna("")
            converted_cols.append((col, "datetime/timedelta→string"))
            continue

        # 3) Pure numeric and boolean are fine
        if pd.api.types.is_bool_dtype(s) or pd.api.types.is_integer_dtype(s) or pd.api.types.is_float_dtype(s):
            # Leave as is
            continue

        # 4) String dtype
        if pd.api.types.is_string_dtype(s):
            # Ensure no non-scalars slipped in
            if s.map(lambda x: isinstance(x, SAFE_SCALARS) or x is None or (isinstance(x, float) and np.isnan(x))).all():
                continue
            # Mixed values in a string dtype series → coerce elementwise
            obs[col] = s.map(_stringify_complex)
            converted_cols.append((col, "string-mixed→stringified"))
            continue

        # 5) Object or odd dtypes
        if pd.api.types.is_object_dtype(s):
            # If every non-null is a safe scalar -> cast to string
            non_null = s.dropna()
            if non_null.map(lambda x: isinstance(x, SAFE_SCALARS)).all():
                obs[col] = s.astype("string").fillna("")
                converted_cols.append((col, "object-safe→string"))
            else:
                # Elementwise stringify
                obs[col] = s.map(_stringify_complex)
                converted_cols.append((col, "object-complex→json-string"))
            continue

    if verbose and converted_cols:
        print("[sanitize_obs_for_h5ad] Converted columns:")
        for name, how in converted_cols:
            print(f"  - {name}: {how}")

    adata.obs = obs
    return adata


def read_features_1or2col(path: str) -> pd.DataFrame:
    """Accepts 1-col or 2-col features tsv.gz. Returns DataFrame with columns ['feature_id','feature_name']."""
    df = pd.read_csv(path, sep="\t", header=None)
    if df.shape[1] == 1:
        df.columns = ["feature_id"]
        df["feature_name"] = df["feature_id"]
    elif df.shape[1] >= 2:
        df = df.iloc[:, :2]
        df.columns = ["feature_id", "feature_name"]
    else:
        raise ValueError(f"Empty features file: {path}")
    df["feature_id"] = df["feature_id"].astype(str)
    df["feature_name"] = df["feature_name"].astype(str)
    return df

def read_barcodes_1col(path: str) -> np.ndarray:
    bc = pd.read_csv(path, sep="\t", header=None)
    if bc.shape[1] != 1:
        raise ValueError(f"Barcodes file must be 1 column: {path}")
    return bc.iloc[:,0].astype(str).values

def read_mtx_triplet(stem: str) -> ad.AnnData:
    """Load MatrixMarket + barcodes + features into AnnData with raw counts. Expects genes x cells in MTX."""
    mtx = io.mmread(stem + ".mtx.gz").tocsr().astype(np.float32)
    barcodes = read_barcodes_1col(stem + "_barcodes.tsv.gz")
    features = read_features_1or2col(stem + "_features.tsv.gz")

    if mtx.shape[0] != features.shape[0]:
        raise ValueError(f"MTX rows {mtx.shape[0]} != features {features.shape[0]} for {stem}")
    if mtx.shape[1] != barcodes.shape[0]:
        raise ValueError(f"MTX cols {mtx.shape[1]} != barcodes {barcodes.shape[0]} for {stem}")

    adata = ad.AnnData(X=mtx.T)  # cells x genes
    adata.obs_names = pd.Index(barcodes, name="barcode").astype(str)
    # enforce unique names for scvi
    adata.obs_names = adata.obs_names.astype(str).map(str)
    adata.var_names = pd.Index(features["feature_id"], name="feature_id").astype(str)
    adata.var_names_make_unique()
    adata.var["feature_name"] = features["feature_name"].values
    return adata

def safe_read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()



def get_cell2loc_abundance_and_long_proportions(adata):
    """
    Returns:
        abund_df : DataFrame of predicted cell abundances (means_cell_abundance_w_sf)
        prop_long_df : long-format DataFrame with columns [spot, cell_type, proportion]
    """
    # 1. Extract predicted abundances
    abund_df = adata.obsm["means_cell_abundance_w_sf"].copy()

    # 2. Clean prefixes
    abund_df.columns = [
        c.replace("meanscell_abundance_w_sf_", "")
         .replace("means_cell_abundance_w_sf_", "")
         .replace("cell_abundance_w_sf_", "")
         .replace("cell_abundance_", "")
        for c in abund_df.columns
    ]

    # 3. Ensure numeric + nonnegative
    abund_df = abund_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    abund_df[abund_df < 0] = 0.0

    # 4. Normalize per spot to get proportions
    prop_df = abund_df.div(abund_df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)

    # 5. Convert to long / tidy format
    prop_long_df = prop_df.reset_index().melt(
        id_vars=prop_df.index.name or "index",
        var_name="cell_type",
        value_name="proportion"
    ).rename(columns={prop_df.index.name or "index": "spot"})

    return abund_df, prop_long_df



# ----------------- Load inputs -----------------

# Load parts
adata_sc = read_mtx_triplet(os.path.join(comp_dir, f"sc_count"))
adata_sp = read_mtx_triplet(os.path.join(comp_dir, "spatial_count"))

# Convert sparse matrices to CSR for scvi/cell2location
if sp.issparse(adata_sc.X): adata_sc.X = adata_sc.X.tocsr()
if sp.issparse(adata_sp.X): adata_sp.X = adata_sp.X.tocsr()

sc_meta = safe_read_csv(os.path.join(comp_dir, "sc_meta.csv"))
sp_meta = pd.read_csv("/project/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Annotations/adata_metadata.csv")
sp_loc = safe_read_csv(os.path.join(comp_dir, "spatial_location.csv"))

# Join spatial metadata
obs_df = adata_sp.obs.reset_index().rename(columns={"index": "barcode"})
if not sp_meta.empty:
    sp_meta = sp_meta.rename(columns={"spot": "barcode"})
    obs_df = obs_df.merge(sp_meta, how="left", on="barcode")
adata_sp.obs = obs_df.set_index("barcode")


# Join sc meta
if not sc_meta.empty:
    if "barcode" in sc_meta.columns:
        sc_meta = sc_meta.set_index("barcode")
        adata_sc.obs = adata_sc.obs.join(sc_meta, how="left")
    else:
        sc_meta = sc_meta.copy()
        sc_meta["barcode"] = adata_sc.obs.index.astype(str)
        sc_meta = sc_meta.set_index("barcode")
        adata_sc.obs = adata_sc.obs.join(sc_meta, how="left")

# Force cell type present
if "CellType_Deconvolution" not in adata_sc.obs.columns:
    raise ValueError("sc_meta.csv must contain 'CellType_Deconvolution'")
adata_sc.obs["CellType_Deconvolution"] = adata_sc.obs["CellType_Deconvolution"].astype("category")

# Force donor present
if "donor_id" not in adata_sc.obs.columns:
    adata_sc.obs["donor_id"] = "donor_1"

# drop unknown / other cell types and enforce >= 25 cells per type
obs = adata_sc.obs.copy()
mask_valid = (
    obs["CellType_Deconvolution"].notna()
    & ~obs["CellType_Deconvolution"].isin(["unknown", "Unknown", "Other", "other"])
)
obs = obs[mask_valid]
counts = adata_sc.X[mask_valid.values, :]

ct_counts = obs["CellType_Deconvolution"].value_counts()
keep_ct = ct_counts[ct_counts >= 25].index
mask_ct = obs["CellType_Deconvolution"].isin(keep_ct)
obs = obs[mask_ct]
counts = counts[mask_ct.values, :]

adata_sc = ad.AnnData(X=counts, obs=obs, var=adata_sc.var.copy())


# Set spatial coordinates
sp_loc["barcode"] = adata_sp.obs.index
sp_loc = sp_loc.set_index("barcode")
adata_sp.obs = adata_sp.obs.join(sp_loc, how="left")


# Check/load RegressionModel 
reg_model_pt = os.path.join(out_dir, "regression_model", "model.pt")
ref_h5ad     = os.path.join(out_dir, "Ref_sc.h5ad")
ref_csv      = os.path.join(out_dir, f"{organ}_Reference_InferredAverages_df.csv")

if os.path.isfile(reg_model_pt) and os.path.exists(ref_h5ad) and os.path.exists(ref_csv):
    print(f"[INFO] Found RegressionModel at: {reg_model_pt}")

    adata_sc = ad.read_h5ad(ref_h5ad)
    
    inf_aver = pd.read_csv(ref_csv, index_col=0)
    inf_aver.index = inf_aver.index.astype(str)

else:
    print(f"[INFO] No RegressionModel found at: {reg_model_pt}. Training now.")

    # ----------------- Filter counts -----------------
    # Drop all-zero genes only
    adata_sc = adata_sc[:, (adata_sc.X.sum(axis=0) > 0).A1].copy()
    adata_sp = adata_sp[:, (adata_sp.X.sum(axis=0) > 0).A1].copy()

    # remove MT genes from var_names
    if "feature_name" not in adata_sp.var.columns:
        raise KeyError("adata_sp.var['feature_name'] is required for MT gene detection.")
    mt_mask = adata_sp.var["feature_name"].str.upper().str.startswith("MT-").fillna(False).values
    X_mt = adata_sp[:, mt_mask].X
    adata_sp.obsm["MT"] = X_mt.toarray() if sp.issparse(X_mt) else np.asarray(X_mt)
    adata_sp = adata_sp[:, ~mt_mask].copy()

    # Filter ref genes
    selected = filter_genes(adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata_sc = adata_sc[:, selected].copy()


    # ----------------- Train RegressionModel on scRNA (raw counts) -----------------
    print(f"[INFO] {organ}: training RegressionModel (scRNA reference)")
    RegressionModel.setup_anndata(
        adata_sc,
        batch_key="donor_id",                  # handle donor/batch effects
        labels_key="CellType_Deconvolution"    # cell type labels in scRNA
    )
    reg = RegressionModel(adata_sc)
    reg.train(max_epochs=250, batch_size=1024)

    adata_sc = reg.export_posterior(
        adata_sc,
        sample_kwargs={"num_samples": 1000, 'batch_size': 2500}
    )

    reg.save(os.path.join(out_dir, "regression_model"), overwrite=True)

    # Save anndata object with results
    adata_sc = sanitize_obs_for_h5ad(adata_sc, verbose=True)
    adata_file = f"{out_dir}/Ref_sc.h5ad"

    for col in adata_sc.obs.columns:
        if pd.api.types.is_string_dtype(adata_sc.obs[col].dtype):
            adata_sc.obs[col] = adata_sc.obs[col].astype(object).astype(str)

    adata_sc.write(adata_file)

    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
        inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}' 
                                        for i in adata_sc.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_sc.uns['mod']['factor_names']
    inf_aver.iloc[0:5, 0:5]

    inf_aver.to_csv(f"{out_dir}/{organ}_Reference_InferredAverages_df.csv", index=True)



# SAMPLE LOOP -----------------------------------------------------------------------------------------------------------


adata_file = f"{out_dir}/{organ}_{sample}_results.h5ad"
out_csv = os.path.join(out_dir, f"{organ}_{sample}_CellAbundance.csv")

# Skip if results already exist
if os.path.exists(adata_file) and os.path.exists(out_csv):
    print(f"Skipping {sample} — already completed.")
    sys.exit(f"[ALREADY COMPLETED]")

if os.path.exists(adata_file) and not os.path.exists(out_csv):
    adata_sp_sub = ad.read_h5ad(adata_file)
    matrix_df, cell_abundance = get_cell2loc_abundance_and_long_proportions(adata_sp_sub)


    # Write to CSV
    out_csv = os.path.join(out_dir, f"{organ}_{sample}_CellAbundance.csv")
    print(f"saving ST proportions csv to {out_csv}")
    cell_abundance.to_csv(out_csv, index=True)


    out_mat = os.path.join(out_dir, f"{organ}_{sample}_CellAbundance_matrix.csv")
    print(f"saving ST proportions matrix csv to {out_mat}")
    matrix_df.to_csv(out_mat, index=True)

    # Verify saved files exist
    for f in [adata_file, out_csv]:
        if not os.path.isfile(f):
            sys.exit(f"[ERROR] File not found after save attempt: {f}")

else:

    adata_sp_sub = adata_sp[adata_sp.obs["Sample"] == sample].copy()

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_sp_sub.var_names.astype(str), inf_aver.index.astype(str))

    if intersect.size == 0:
        raise ValueError(f"No shared genes between spatial {sample} and reference signatures.")

    intersect = pd.Index(intersect).astype(str)
    adata_sp_sub = adata_sp_sub[:, intersect].copy()
    inf_aver_sample = inf_aver.loc[intersect, :].copy()

    # prepare cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_sp_sub)

    mod = cell2location.models.Cell2location(
        adata_sp_sub, 
        cell_state_df=inf_aver_sample, 
        N_cells_per_location=10,
        detection_alpha=20
    ) 
    scvi.settings.dl_num_workers = 1  # avoids NFS thrash on Midway


    mod.train(
        max_epochs=30000,
        batch_size=None,
        train_size=1.0
    )


    # export the estimated cell abundance (summary of the posterior distribution).
    adata_sp_sub = mod.export_posterior(
        adata_sp_sub,
        sample_kwargs={"num_samples": 300, "batch_size": 2000}
    )
    print(f"saving mod to {out_dir}/{organ}_{sample}")
    mod.save(f"{out_dir}/{organ}_{sample}", overwrite=True)

    # Save anndata object with results

    if os.path.exists(adata_file) and os.path.exists(out_csv):
        print(f"Skipping {sample} — already completed.")
        sys.exit(f"[ALREADY COMPLETED]")


    print(f"saving ST adata to {out_dir}/{organ}_{sample}_results.h5ad")
    adata_sp_sub = sanitize_obs_for_h5ad(adata_sp_sub, verbose=True)
    adata_file = f"{out_dir}/{organ}_{sample}_results.h5ad"

    for col in adata_sp_sub.obs.columns:
        if pd.api.types.is_string_dtype(adata_sp_sub.obs[col].dtype):
            adata_sp_sub.obs[col] = adata_sp_sub.obs[col].astype(object).astype(str)

    adata_sp_sub.write(adata_file)

    matrix_df, cell_abundance = get_cell2loc_abundance_and_long_proportions(adata_sp_sub)



    # Write to CSV
    out_csv = os.path.join(out_dir, f"{organ}_{sample}_CellAbundance.csv")
    print(f"saving ST proportions csv to {out_csv}")
    cell_abundance.to_csv(out_csv, index=True)


    out_mat = os.path.join(out_dir, f"{organ}_{sample}_CellAbundance_matrix.csv")
    print(f"saving ST proportions matrix csv to {out_mat}")
    matrix_df.to_csv(out_mat, index=True)

    # Verify saved files exist
    for f in [adata_file, out_csv]:
        if not os.path.isfile(f):
            sys.exit(f"[ERROR] File not found after save attempt: {f}")

    

