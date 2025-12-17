# Setup -----------------------------------------------------------

source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")



# FUNCTIONS ----------------

clean_text = function(x) {
  toupper(gsub(" ", "", gsub("[^A-Za-z0-9 ]", "", x)))
}

missing_broad = function(df){
  if (any(is.na(df$CellType_Broad_Deconvolution))) {
    missing_rows = which(is.na(df$CellType_Broad_Deconvolution))
    stop(
      paste0(
        "Need to harmonize! Missing values in CellType_Broad_Deconvolution at rows: ",
        paste(missing_rows, collapse = ", "),
        "\nBarcodes: ",
        paste(df$CellType_Deconvolution[missing_rows], collapse = ", ")
      )
    )
  }
} 

ref_meta_read = function(ref_dir, organ){
    read.csv2(file.path(ref_dir, organ, "components", "sc_meta.csv")) %>%
    mutate(
      CellType_Broad_Deconvolution = case_when(
        clean_text(CellType_Deconvolution) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType_Deconvolution)],
        TRUE ~ NA
      )
    )
}

get_top5_long = function(pred_df,
                          id_cols = c("barcode","y","x","Sample","Organ"),
                          k = 5) {
  stopifnot(all(id_cols %in% names(pred_df)))
  stopifnot(k >= 1)
  
  prob_cols = setdiff(names(pred_df), id_cols)
  stopifnot(length(prob_cols) > 0)
  
  pred_df %>%
    mutate(.row_id = row_number()) %>%
    pivot_longer(
      cols = all_of(prob_cols),
      names_to = "CellType_Deconvolution",
      values_to = "proportion"
    ) %>%
    group_by(.row_id, across(all_of(id_cols))) %>%
    arrange(desc(proportion), CellType_Deconvolution, .by_group = TRUE) %>%
    mutate(.rank = row_number()) %>%
    filter(.rank <= k) %>%
    ungroup() %>%
    select(all_of(id_cols), CellType_Deconvolution, proportion)
}



# Load REFERENCE DATA --------------

Metadata = read.csv2("Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/IHC/dfs/adata_obs.csv")

Comp_dictionary = read.csv2("Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Dictionaries/Broad_dictionary.csv") %>%
  deframe

names(Comp_dictionary) = clean_text(names(Comp_dictionary))

ref_dir = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Reference_data"

organs = c("LU", "SP", "LI", "KI")



# Compile CARD

CARD_files = list.files(
  "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/CARD_output",
  full.names = T,
  recursive = T,
  pattern = "CARD_proportions\\.csv$"
)


CARD_files_df = map_dfr(CARD_files, function(file){
  data.frame(
    Organ = basename(dirname(file)),
    path = file
  )
}) %>%
  arrange(Organ, path)

CARD_df = map_dfr(CARD_files_df$path, function(file){
  
  df = read.csv2(file)
  
  organ = unique(df$Organ)
  
  ref_meta = ref_meta_read(ref_dir, organ)
    
  missing_broad(ref_meta)
  
  possible_celltype_broad = unique(ref_meta$CellType_Broad_Deconvolution)
  
  possible_spots = Metadata  %>%
    filter(Organ == organ, CellType_Broad %in% possible_celltype_broad)
  
  df$spot = NULL
  
  df_return = df %>% 
    dplyr::rename(spot = barcode) %>%
    filter(spot %in% possible_spots$spot) %>%
    group_by(spot) %>%
    arrange(desc(proportion)) %>%
    dplyr::slice(1:3) %>%
    ungroup() %>%
    mutate(
      CellType_Broad_Deconvolution = case_when(
        clean_text(CellType_Deconvolution) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType_Deconvolution)],
        TRUE ~ NA
      ),
      Decon_Method = "Card"
    )  %>%
    select(spot, x, y, CellType_Deconvolution, CellType_Broad_Deconvolution, proportion, Decon_Method) %>%
    left_join(
      .,
      Metadata %>% 
        select(spot, Organ, Subregion, CellType, CellType_Broad, Sample, Treatment),
      by = "spot"
    ) %>%
    mutate(
      Match = case_when(
        CellType_Broad_Deconvolution == CellType_Broad ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    group_by(spot) %>%
    arrange(desc(Match == "Yes"), desc(proportion), .by_group = TRUE) %>% 
    dplyr::slice(1) %>%
    ungroup() 
  
  missing_broad(df_return)
  
  df_return
  
})


# Compile RCTD ------------------------------------------------------------------------------------------------------------------------------



RCTD_files = list.files(
  "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/RCTD_output",
  full.names = T,
  recursive = T,
  pattern = "RCTD_proportions\\.csv$"
)


RCTD_files_df = map_dfr(RCTD_files, function(file){
  data.frame(
    Organ = basename(dirname(file)),
    path = file
  )
}) %>%
  arrange(Organ, path)


RCTD_df = map_dfr(RCTD_files_df$path, function(file){
  
  df = read.csv2(file)
  
  if("CellType_Deconvolution" %!in% colnames(df)){
    
    if("Sample.Sample" %in% colnames(df)){
      df = df %>% dplyr::rename(Sample = 'Sample.Sample')
    }
    if("Organ.Organ" %in% colnames(df)){
      df = df %>% dplyr::rename(Organ = 'Organ.Organ')
    }
    
    df = get_top5_long(df)
  }
  
  organ = unique(df$Organ)
  
  
  ref_meta = ref_meta_read(ref_dir, organ)
  
  
  missing_broad(ref_meta)
  
  possible_celltype_broad = unique(ref_meta$CellType_Broad_Deconvolution)
  
  possible_spots = Metadata  %>%
    filter(Organ == organ, CellType_Broad %in% possible_celltype_broad)
  
  df$spot = NULL
  df_return = df %>%
    filter(!is.na(CellType_Deconvolution)) %>%
    dplyr::rename(spot = barcode) %>%
    filter(spot %in% possible_spots$spot) %>%
    group_by(spot) %>%
    arrange(desc(proportion)) %>%
    dplyr::slice(1:3) %>%
    ungroup() %>%
    mutate(
      CellType_Broad_Deconvolution = case_when(
        clean_text(CellType_Deconvolution) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType_Deconvolution)],
        TRUE ~ NA
      ),
      Decon_Method = "RCTD"
    ) %>%
    select(spot, x, y, CellType_Deconvolution, CellType_Broad_Deconvolution, proportion, Decon_Method) %>%
    left_join(
      .,
      Metadata %>% 
        select(spot, Organ, Subregion, CellType, Sample, Treatment) %>%
        mutate(
          CellType_Broad = case_when(
            clean_text(CellType) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType)],
            TRUE ~ NA
          )
        ),
      by = "spot"
    ) %>%
    mutate(
      Match = case_when(
        CellType_Broad_Deconvolution == CellType_Broad ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    group_by(spot) %>%
    arrange(desc(Match == "Yes"), desc(proportion), .by_group = TRUE) %>% 
    dplyr::slice(1) %>%
    ungroup() 
  
  missing_broad(df_return)
  
  
  df_return
  
})





# Compile cell2location --------------------------------------------------------------------------------------------------------------------


cell2location_files = list.files(
  "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/cell2location_output",
  full.names = T,
  recursive = T,
  pattern = "_[12]_CellAbundance\\.csv$"
)


cell2location_files_df = map_dfr(cell2location_files, function(file){
  data.frame(
    Organ = basename(dirname(file)),
    Sample = toupper(regmatches(basename(file),
                                regexpr("(CTRL|LPS)_\\d+", basename(file), ignore.case = TRUE))),
    path = file
  )
}) %>%
  arrange(Organ, Sample, path) 

Cell2Loc_df = map_dfr(cell2location_files_df$path, function(file){
  
  organ = dirname(file) %>% basename
  df = read.csv2(file)
  
  ref_meta = ref_meta_read(ref_dir, organ)
  
  missing_broad(ref_meta)
  
  
  possible_celltype_broad = unique(ref_meta$CellType_Broad_Deconvolution)
  
  possible_spots = Metadata  %>%
    filter(Organ == organ, CellType_Broad %in% possible_celltype_broad)
  
  if("X" %in% colnames(df)){ df$X = NULL }
  df_return = df %>%
    dplyr::rename(CellType_Deconvolution = cell_type) %>%
    filter(!is.na(CellType_Deconvolution)) %>%
    filter(spot %in% possible_spots$spot) %>%
    group_by(spot) %>%
    arrange(desc(proportion)) %>%
    dplyr::slice(1:3) %>%
    ungroup() %>%
    mutate(
      CellType_Broad_Deconvolution = case_when(
        clean_text(CellType_Deconvolution) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType_Deconvolution)],
        TRUE ~ NA
      ),
      Decon_Method = "Cell2Location"
    ) %>%
    select(spot, CellType_Deconvolution, CellType_Broad_Deconvolution, proportion, Decon_Method) %>%
    left_join(
      .,
      Metadata %>% 
        dplyr::rename(x = x_scaled_image, y = y_scaled_image) %>%
        select(spot, Organ, Subregion, CellType, Sample, Treatment, x, y) %>%
        mutate(
          CellType_Broad = case_when(
            clean_text(CellType) %in% clean_text(names(Comp_dictionary)) ~ Comp_dictionary[clean_text(CellType)],
            TRUE ~ NA
          )
        ),
      by = "spot"
    ) %>%
    mutate(
      Match = case_when(
        CellType_Broad_Deconvolution == CellType_Broad ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    group_by(spot) %>%
    arrange(desc(Match == "Yes"), desc(proportion), .by_group = TRUE) %>% 
    dplyr::slice(1) %>%
    ungroup() %>%
    select(
      spot, x, y, CellType_Deconvolution, CellType_Broad_Deconvolution, proportion, Decon_Method, Organ,                       
      Subregion, CellType_Broad, CellType, Sample, Treatment, Match  
    )
  
  missing_broad(df_return)
  
  
  df_return
  
})








# Compile ---------------------------------------------------------------------------------------------------------------------------------

full_df = rbind(
  RCTD_df,
  CARD_df,
  Cell2Loc_df
)


CellType_Deconvolution_df = rbind(
  Metadata %>% 
    dplyr::rename(x = x_scaled_image, y = y_scaled_image) %>%
    mutate(
      CellType_Deconvolution = CellType,
      CellType_Deconvolution_Broad = CellType_Broad,
      Decon_Method = "CellKB",
      proportion = 1,
      Match = "Yes"
    ) %>%
    select(all_of(colnames(full_df))) %>%
    filter(spot %in% full_df$spot),
  full_df
)

write.csv(
  CellType_Deconvolution_df,
  file = "Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Annotations/deconvolution_results.csv",
  row.names = F
)

