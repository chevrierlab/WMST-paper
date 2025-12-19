# SETUP -------------------

source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")

# ---------------------------------


# PERCENTER ---------------------------------------------------------------------

percenter = function(x,y){round(as.numeric(x)/as.numeric(y)*100,0)}



# COLOR PALETTES ------------------------------------------------------

core_celltypes = unique(regulons_df$CellType_Core) %>% sort()

m_color_palette = c(
  "#B76674", 
  "#B4866A", 
  "#9C9A4A", 
  "#5F8271",  
  "#5D84A2", 
  "#6D5691", 
  "#A36678",  
  "#5E7C5B",  
  "#3F6194"  
)

core_palette = setNames(m_color_palette, core_celltypes)


generate_greys =   function(n, start_hex = "#EAF2FB", end_hex = "#394E6A", space = "Lab") {
  stopifnot(n > 0)
  pal  = grDevices::colorRampPalette(c(start_hex, end_hex), space = space)
  cols = pal(n)
  names(cols) = as.character(seq_len(n))
  cols
}


# SIZING ----------------------------------------------------------

title_size = 8
legend_size = 7
column_width = 4.8
column_width_secondary = 2.8


# GENERATE AND SAVE HM -----------------------


Save_Final_Heatmaps = function(HM1, HM2, HM3, title, out_dir, h, w){
  
  HM_Membership = HM1
  HM_CoreDE_BW = HM2
  HM_CoreDE_RB = HM3
  
  cluster_df = HM_Membership$cluster_df
  
  heatmaps = list(
    Membership = HM_Membership$heatmap,
    BW         = HM_CoreDE_BW$heatmap,
    RB         = HM_CoreDE_RB$heatmap
  )
  
  combos = list(
    "full_HM"              = c("Membership", "BW", "RB")
  )
  
  save_heatmap_combos(
    title = paste0(title),
    out_dir = folder.path(out_dir),
    heatmaps = heatmaps,
    combos = combos,
    width_in = w,
    height_in = h,
    res = 900,
    n_genes = length(unique(cluster_df$gene))
  )
  
  
  # Save specs
  clean_title = str_replace_all(title, " ", "_") %>%
    str_replace_all(., ",", "") %>%
    str_replace_all(., ":", "")
  
  saveRDS(HM_Membership$matrix, file = folder.path(out_dir, paste0(clean_title, "_HM_Info"), "Membership_matrix.rds"), compress = "xz")
  saveRDS(HM_CoreDE_BW$matrix, file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "BinaryDE_matrix.rds"), compress = "xz")
  saveRDS(HM_CoreDE_RB$matrix, file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "LFC2_matrix.rds"), compress = "xz")
  
  saveRDS(HM_Membership$heatmap, file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "Membership_hm.rds"), compress = "xz")
  saveRDS(HM_CoreDE_BW$heatmap, file = folder.path(out_dir, paste0(clean_title, "_HM_Info"), "BinaryDE_hm.rds"), compress = "xz")
  saveRDS(HM_CoreDE_RB$heatmap, file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "LFC2_hm.rds"), compress = "xz")
  
  write.csv(
    cluster_df,
    file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "Cluster_df.csv"),
    row.names = F
  )
  
  write.csv(
    as.tibble(HM_Membership$cols),
    file = folder.path(out_dir,  paste0(clean_title, "_HM_Info"), "Column_order.csv"),
    row.names = F
  )
  
}


# SAVE FUNCTIONS --------------------------------------------------


save_heatmap_combos = function(
    title,
    n_genes,
    out_dir,
    heatmaps,
    combos,
    height_in = 6.71,
    width_in = 5.67,
    res = 600,
    merge_legend = TRUE
) {
  stopifnot(dir.exists(out_dir))
  stopifnot(is.list(heatmaps), length(heatmaps) > 0)
  stopifnot(is.list(combos), length(combos) > 0)
  
  assemble_plus = function(keys) {
    vals = lapply(keys, function(k){
      if (!k %in% names(heatmaps)) stop(sprintf("Heatmap key '%s' missing.", k))
      heatmaps[[k]]
    })
    Reduce(`+`, vals)
  }
  
  for (nm in names(combos)) {
    hm_keys = combos[[nm]]
    combo = assemble_plus(hm_keys)
    

    clean_title = str_replace_all(title, " ", "_")%>%
      str_replace_all(., ",", "") %>%
      str_replace_all(., ":", "")
    
    
    save_heatmap(
      file_path = file.path(out_dir, paste0(clean_title, "_", nm, ".png")),
      draw_code = {
        draw(
          combo,
          merge_legend = merge_legend,
          heatmap_legend_side = "right",
          column_title = paste0(title, " (", n_genes, ")"),
          column_title_side = "top",
          column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
          padding = grid::unit(c(2, 2, 8, 2), "mm")
        )
      },
      width_in = width_in,
      height_in = height_in,
      res = res
    )
    
    save_heatmap(
      file_path = file.path(out_dir, paste0(clean_title, "_", nm, ".pdf")),
      draw_code = {
        draw(
          combo,
          merge_legend = merge_legend,
          heatmap_legend_side = "right",
          column_title = paste0(title, " (", n_genes, ")"),
          column_title_side = "top",
          column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
          padding = grid::unit(c(2, 2, 8, 2), "mm")
        )
      },
      width_in = width_in,
      height_in = height_in,
      res = res
    )
  }
}



ht_isolate = function(expr) {
  old = ht_opt(RESET = TRUE)
  on.exit(ht_opt(old), add = TRUE)
  ht_opt(save_last = FALSE, verbose = FALSE)
  
  try(grid::seekViewport("ROOT"), silent = TRUE)
  grid::upViewport(0)
  grid::grid.newpage()
  
  force(expr)
}

close_all_devices = function() {
  dl = dev.list()
  if (!is.null(dl)) {
    while (!is.null(dev.list())) try(dev.off(), silent = TRUE)
  }
  invisible(NULL)
}

save_heatmap = function(file_path, draw_code,
                        width_in = 5.67, height_in = 8.69,
                        res = 300) {
  
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  
  ext = tolower(tools::file_ext(file_path))
  if (ext == "png") {
    grDevices::png(filename = file_path, width = width_in, height = height_in,
                   units = "in", res = res, type = "cairo")
  } else if (ext == "pdf") {
    grDevices::pdf(file = file_path, width = width_in, height = height_in, onefile = TRUE)
  } else if (ext == "svg") {
    grDevices::svg(file = file_path, width = width_in, height = height_in)
  } else {
    stop("Unsupported extension: ", ext)
  }
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  grid::grid.newpage()
  
  old = ComplexHeatmap::ht_opt(RESET = TRUE)
  on.exit(ComplexHeatmap::ht_opt(old), add = TRUE)
  ComplexHeatmap::ht_opt(save_last = FALSE, verbose = FALSE)
  
  force(draw_code)
  
  invisible(TRUE)
}




# PATHWAY_ENRICH_UNIVERSE -------------------------------------------------------------

safe_enrich = function(genes, db,
                       tries = 14L,
                       pause = 0.8,
                       bases = c("https://maayanlab.cloud/Enrichr",
                                 "https://amp.pharm.mssm.edu/Enrichr")) {
  genes = unique(stats::na.omit(genes))
  
  empty_tbl = function() {
    tibble::tibble(
      Term = character(),
      Overlap = character(),
      Adjusted.P.value = numeric(),
      Combined.Score = numeric(),
      Genes = character()
    )
  }
  
  if (length(genes) == 0) return(empty_tbl())
  
  for (i in seq_len(tries)) {
    
    if (length(bases) > 0 && requireNamespace("enrichR", quietly = TRUE)) {
      enrichR::setEnrichrSite(bases[(i - 1) %% length(bases) + 1])
    }
    
    res = try(enrichr(genes, databases = db)[1], silent = TRUE)
    
    if (!inherits(res, "try-error") && is.list(res) && !is.null(res[[db]])) {
      out = tibble::as_tibble(res[[db]])
      

      if ("genes" %in% names(out) && !"Genes" %in% names(out)) out = dplyr::rename(out, Genes = genes)
      if ("overlap" %in% names(out) && !"Overlap" %in% names(out)) out = dplyr::rename(out, Overlap = overlap)
      
      need = c("Term","Overlap","Adjusted.P.value","Combined.Score","Genes")
      missing = setdiff(need, names(out))
      if (length(missing) == 0) return(dplyr::select(out, dplyr::all_of(need)))
      
      return(empty_tbl())
    }
    
    Sys.sleep(pause * 2^(i - 1))
  }
  
  empty_tbl()
}



safe_enrich_cp = function(genes,
                          universe = NULL,
                          p_cutoff = 0.05,
                          q_cutoff = 0.05,
                          min_gs = 10,
                          max_gs = 500) {
  stopifnot(is.character(genes))
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("Install clusterProfiler and org.Mm.eg.db")
  }
  
  tryCatch({
    eg = clusterProfiler::enrichGO(
      gene           = unique(genes),
      universe       = if (is.null(universe)) NULL else unique(universe),
      OrgDb          = org.Mm.eg.db::org.Mm.eg.db,
      keyType        = "SYMBOL",
      ont            = "BP",
      pAdjustMethod  = "BH",
      pvalueCutoff   = p_cutoff,
      qvalueCutoff   = q_cutoff,
      minGSSize      = min_gs,
      maxGSSize      = max_gs,
      readable       = FALSE
    )
    
    if (is.null(eg) || nrow(as.data.frame(eg@result)) == 0) return(tibble::tibble())
    
    
    out = dplyr::as_tibble(eg@result) |>
      dplyr::mutate(
        GeneRatioNum = vapply(GeneRatio, function(x){p = strsplit(x, "/", fixed = TRUE)[[1]]; as.numeric(p[1])/as.numeric(p[2])}, numeric(1)),
        BgRatioNum   = vapply(BgRatio,   function(x){p = strsplit(x, "/", fixed = TRUE)[[1]]; as.numeric(p[1])/as.numeric(p[2])}, numeric(1))
      ) |>
      dplyr::arrange(p.adjust)
    
    print(nrow(out))
    
    out = out %>%
      dplyr::rename(
        Term = Description,
        Overlap = GeneRatio,
        P.value = pvalue,
        Adjusted.P.value = p.adjust,
        Genes = geneID
      ) %>%
      mutate(DB = "GO") %>%
      select(ID, Term, Overlap, P.value, Adjusted.P.value, Genes)
    
    print(nrow(out))
    return(out)
  }, error = function(e) tibble::tibble())
}


pathway_enrich_UNIVERSE = function(genes, universe, cluster){
  

  output_GO = safe_enrich_cp(genes, universe = universe)
  output_wiki = safe_enrich(genes, "WikiPathways_2024_Mouse")
  
  
  top_n = 30  
  
  GO_top = output_GO %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::mutate(DB = "GO_Biological_Process_2025") %>%
    dplyr::mutate(
      rows = ID
    ) %>%
    tibble::column_to_rownames("rows")
  
  Wiki_top = output_wiki %>%
    dplyr::arrange(Adjusted.P.value) %>%
    dplyr::slice(1:top_n) %>%
    dplyr::mutate(
      ID   = stringr::str_extract(Term, "WP\\d+"),
      Term = stringr::str_trim(stringr::str_remove(Term, "\\s+WP\\d+$")),
      rows = ID
    ) %>%
    tibble::column_to_rownames("rows") %>%
    dplyr::mutate(DB = "WikiPathways_2024_Mouse")
  
  Combined_top = rbind(
    Wiki_top %>% dplyr::select(DB, ID, Term, Overlap, Adjusted.P.value, Genes),
    GO_top   %>% dplyr::select(DB, ID, Term, Overlap, Adjusted.P.value, Genes)
  )
  
  Out_GO = GO_top[1:4, ] %>%
    dplyr::mutate(source = "GO_BP") %>%
    dplyr::select(source, Term, Adjusted.P.value)
  
  Out_Wiki = Wiki_top[1:3, ] %>%
    dplyr::mutate(source = "Wiki") %>%
    dplyr::select(source, Term, Adjusted.P.value)
  
  cluster_identifier = data.frame(
    source = "identifier",
    Term = paste0("CLUSTER ", cluster, ":     "),
    Adjusted.P.value=0
  )
  
  complete_Out = rbind(cluster_identifier, Out_GO, Out_Wiki)
  
  
  Out_GO2 = GO_top[1:10, ] %>%
    dplyr::mutate(source = "GO_BP") %>%
    dplyr::select(source, Term, Adjusted.P.value)
  
  Out_Wiki2 = Wiki_top[1:10, ] %>%
    dplyr::mutate(source = "Wiki") %>%
    dplyr::select(source, Term, Adjusted.P.value)
  
  
  complete_Out2 = rbind(Out_Wiki2, Out_GO2)
  
  
  list(
    df       = Combined_top,
    selected = complete_Out,
    more_selected = complete_Out2,
    go_bp_out= Out_GO2
  )
}




# CORE_LFC_HM -----------------------------------------------------------------------


Core_LFC_HM = function(core_df, reg_df, row_order, col_order, lfc_max =2){
  
  #  membership per gene × CellType_Core
  membership_df_present =
    core_df %>%
    distinct(CellType_Core, gene, .keep_all = T) %>%
    mutate(value = pmax(-lfc_max, pmin(logfc, lfc_max))) %>%
    select(CellType_Core, gene, value)
  
  membership_df = map_dfr(unique(membership_df_present$CellType_Core), function(cell){
    df = membership_df_present %>%
      filter(CellType_Core == cell)
    df_complete = df %>%
      rbind(
        .,
        reg_df %>%
          dplyr::rename(gene = gene) %>%
          distinct(gene) %>%
          filter(gene %!in% df$gene) %>%
          select(gene) %>%
          mutate(value = 0,
                 CellType_Core = cell) %>%
          select(all_of(colnames(df)))
      ) %>%
      ungroup()
  })
  
  # wide df
  wide_df = 
    membership_df %>%
    pivot_wider(
      id_cols = gene,
      names_from = CellType_Core,
      values_from = value
    )
  
  
  # Convert to matrix
  mat_char =
    wide_df %>%
    select(-gene) %>%
    as.data.frame()
  rownames(mat_char) = wide_df$gene
  mat_char = as.matrix(mat_char)
  col_order = col_order[col_order %in% colnames(mat_char)]
  
  # HM colors
  hm_cols = circlize::colorRamp2(
    c(-lfc_max, 0, lfc_max),
    c("#3b438f", "white", "#c94238") 
  )
  
  
  # Column anno
  
  
  if(ncol(wide_df) ==2){
    final_matrix = mat_char[row_order,col_order] %>%
      as.matrix
    colnames(final_matrix) = colnames(wide_df[2])
  } else{
    final_matrix = mat_char[row_order,col_order]
  }
  
  core_iterations = colnames(final_matrix)
  
  core_factor = factor(core_iterations, levels = unique(core_iterations))
  stopifnot(ncol(final_matrix) == length(core_factor))
  
  core_colors = data.frame(core_palette) %>%
    rownames_to_column("Core_CellTypes") %>%
    filter(Core_CellTypes %in% levels(core_factor)) %>%
    arrange(match(Core_CellTypes, unique(core_iterations))) %>%
    deframe()
  
  
  core_ann = HeatmapAnnotation(
    which = "column",
    blocks = anno_block(
      gp = gpar(fill = core_colors, col = "black"),
    ),
    height = unit(2, "mm")
  )
  
  broad_banner = HeatmapAnnotation(
    which = "column",
    "DE" = rep("all", ncol(final_matrix)),
    col = list("DE" = c(all = "#a5a6d1")),
    show_legend = F,
    simple_anno_size = unit(5, "mm"),
    border = TRUE,
    show_annotation_name = T
  )
  
  # Draw heatmap with split and block annotation
  ht = Heatmap(
    final_matrix,
    name = "LogFC",
    col = hm_cols,
    column_split = core_factor,
    cluster_rows = T,
    cluster_columns = F,
    top_annotation = core_ann,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = T,
    column_names_side = "top",
    column_title = NULL,
    show_column_dend = F,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(-lfc_max, 0, lfc_max),                    
      labels = c(paste0("-",lfc_max), "0", lfc_max),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"),
      labels_gp = gpar(fontsize = legend_size)
    ),
    border = "black",
    width = unit(column_width_secondary , "mm") * ncol(final_matrix)
  )
  
  
  return(
    list(
      heatmap = ht,
      matrix = final_matrix
    )
  )
}




# CORE_DEG_HM -----------------------------------------------------------------------


Core_DEG_HM = function(core_df, reg_df, row_order, col_order){
  
  #  membership per gene × CellType_Core
  membership_df_present =
    core_df %>%
    distinct(CellType_Core, gene) %>%
    mutate(value = 1)
  
  
  membership_df = map_dfr(unique(membership_df_present$CellType_Core), function(cell){
    df = membership_df_present %>%
      filter(CellType_Core == cell)
    df_complete = df %>%
      rbind(
        .,
        reg_df %>%
          dplyr::rename(gene = gene) %>%
          distinct(gene) %>%
          filter(gene %!in% df$gene) %>%
          select(gene) %>%
          mutate(value = 0,
                 CellType_Core = cell) %>%
          select(all_of(colnames(df)))
      ) %>%
      ungroup()
  })
  
  # wide df
  wide_df = 
    membership_df %>%
    pivot_wider(
      id_cols = gene,
      names_from = CellType_Core,
      values_from = value
    )
  
  
  # Convert to matrix
  mat_char =
    wide_df %>%
    select(-gene) %>%
    as.data.frame()
  rownames(mat_char) = wide_df$gene
  mat_char = as.matrix(mat_char)
  col_order = col_order[col_order %in% colnames(mat_char)]
  
  # HM colors
  hm_cols = c("0" = "white", "1" = "black")  
  
  # Column anno
  
  
  if(ncol(wide_df) ==2){
    final_matrix = mat_char[row_order,col_order] %>%
      as.matrix
    colnames(final_matrix) = colnames(wide_df[2])
  } else{
    final_matrix = mat_char[row_order,col_order]
  }
  
  core_iterations = colnames(final_matrix)
  
  core_factor = factor(core_iterations, levels = unique(core_iterations))
  stopifnot(ncol(final_matrix) == length(core_factor))
  
  core_colors = data.frame(core_palette) %>%
    rownames_to_column("Core_CellTypes") %>%
    filter(Core_CellTypes %in% levels(core_factor)) %>%
    arrange(match(Core_CellTypes, unique(core_iterations))) %>%
    deframe()
  
  
  core_ann = HeatmapAnnotation(
    which = "column",
    blocks = anno_block(
      gp = gpar(fill = core_colors, col = "black"),
    ),
    height = unit(2, "mm")
  )
  
  broad_banner = HeatmapAnnotation(
    which = "column",
    "DE" = rep("all", ncol(final_matrix)),
    col = list("DE" = c(all = "#a5a6d1")),
    show_legend = F,
    simple_anno_size = unit(5, "mm"),
    border = TRUE,
    show_annotation_name = T
  )
  
  # 3) Draw heatmap with split and block annotation
  ht = Heatmap(
    final_matrix,
    name = "Significance",
    col = hm_cols,
    column_split = core_factor,
    cluster_rows = F,
    cluster_columns = F,
    top_annotation = core_ann,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = T,
    column_names_side = "top",
    column_title = NULL,
    show_column_dend = F,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(1),
      labels = c("Significant"),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"),   # legend title size
      labels_gp = gpar(fontsize = legend_size)     # legend element label size
    ),
    border = "black",
    width = unit(column_width_secondary , "mm") * ncol(final_matrix)
  )
  
  
  return(
    list(
      heatmap = ht,
      matrix = final_matrix
    )
  )
}




# MEMBERSHIP HEATMAP BUCKETS ---------------------------------

Membership_HM_buckets = function(reg_df, core_df, min_bucket_size = 15){
  
  # membership per gene×CellType_Core 
  membership_df =
    reg_df %>%
    filter(TF %in% c("Stat1", "Irf1")) %>%
    distinct(Organ, CellType_Core, gene = gene, TF) %>%
    mutate(is_Stat1 = TF == "Stat1",
           is_Irf1  = TF == "Irf1") %>%
    group_by(Organ, CellType_Core, gene) %>%
    summarize(has_Stat1 = any(is_Stat1),
              has_Irf1  = any(is_Irf1),
              .groups = "drop") %>%
    mutate(cat = case_when(
      has_Stat1 & has_Irf1 ~ "Both",
      has_Stat1 & !has_Irf1 ~ "Stat1_only",
      !has_Stat1 & has_Irf1 ~ "Irf1_only",
      TRUE ~ "none"
    ))
  
  
  # wide matrix
  wide_df =
    membership_df %>%
    select(Organ, CellType_Core, gene, cat) %>%
    pivot_wider(names_from = CellType_Core, values_from = cat, values_fill = "none")
  
  # Character matrix 
  mat_char =
    wide_df %>%
    select(-Organ, -gene) %>%
    as.data.frame()
  rownames(mat_char) = wide_df$gene
  mat_char = as.matrix(mat_char)
  
  cat_to_num = c(none = 0, Stat1_only = 1, Irf1_only = 2, Both = 3)
  num_mat = matrix(cat_to_num[mat_char],
                   nrow = nrow(mat_char), ncol = ncol(mat_char),
                   dimnames = dimnames(mat_char))
  
  # Heatmap colors
  hm_cols =  c("0" = "#FFFFFF", "1" = "#B86A6A", "2" = "#5A7DA6", "3" = "#6EA27A") # none, Stat1, Irf1, Both
  
  # top annotation
  col_order = colnames(num_mat)
  col_colors = core_palette[names(core_palette) %in% col_order]
  
  
  top_ann = HeatmapAnnotation(
    Organ = anno_block(
      gp = gpar(fill = col_colors, col = "black")
    ),
    height = unit(2, "mm")
  )
  
  # Build clustering
  counts = membership_df %>% dplyr::count(CellType_Core) %>% arrange(desc(n))
  res = build_membership_buckets(membership_df, celltype_priority = counts$CellType_Core, min_n = min_bucket_size)
  Cluster_df1 = res$cluster_df
  
  Cluster_df2 = map_dfr(unique(Cluster_df1$cluster), function(clust){
    
    df_slice = Cluster_df1 %>% filter(cluster == clust)
    
    df_reorder = map_dfr(unique(df_slice$Type), function(type){
      
      genes = df_slice %>% filter(orig_Type == type) %>% pull(gene) %>% as.vector
      
      if(length(genes)==1){
        return( df_slice %>% filter(orig_Type == type) )
      }
      
      mat_slice = num_mat[genes,res$col_order]
      gene_order = HM_orderSlices(mat_slice)
      
      Cluster_df1 %>% filter(cluster == clust) %>% arrange(match(gene, gene_order))
      
    })
    
    
  })
  
  col_order = unique(Cluster_df2$CellType_Core)
  Core_reference = Core_LFC_HM(
    core_df=core_df, 
    reg_df=reg_df, 
    row_order = Cluster_df2$gene, 
    col_order = col_order 
  )
  Core_reference_matrix = Core_reference$matrix
  Membership_reference_matrix = num_mat[Cluster_df2$gene,col_order]
  
  
  
  Cluster_df = map_dfr(unique(Cluster_df2$cluster), function(clust){
    
    df_slice = Cluster_df2 %>% filter(cluster == clust)
    genes = df_slice$gene
    
    mat_slice_df = Membership_reference_matrix[genes,] %>%
      as.df %>%
      mutate(
        summary = do.call(paste, c(.[, col_order, drop = FALSE], sep = "|")),
        Type = consecutive_id(summary)
      )
    
    df_reorder = map_dfr(unique(mat_slice_df$Type), function(type){
      
      genes2 = mat_slice_df %>% filter(Type == type) %>% row.names
      
      if(length(genes)<3){
        return( df_slice %>% filter(gene %in% genes2) )
      }
      
      mat_slice_lfc = Core_reference_matrix[genes2,]
      gene_order = HM_orderSlices(mat_slice_lfc)
      
      df_slice %>% filter(gene %in% genes2) %>% arrange(match(gene, gene_order))
      
    })
    
    
  })
  
  
  
  
  col_split = factor(col_order, levels = col_order)
  
  
  Final_matrix = num_mat[Cluster_df$gene,col_order]
  
  
  cluster_factor = factor(Cluster_df$cluster, levels = unique(Cluster_df$cluster))
  stopifnot(nrow(Final_matrix) == length(cluster_factor))
  
  
  # Cluster anno
  cluster_cols = generate_greys(
    max(as.numeric(Cluster_df$cluster))
  )
  
  cluster_ann = rowAnnotation(
    Cluster = anno_block(
      gp = gpar(fill = unname(cluster_cols), col = "black"),
      labels = levels(cluster_factor),
      labels_gp = gpar(col = "black", fontsize = 12)
    ),
    width = unit(5, "mm")
  )
  
  # Heatmap
  
  ht_final = Heatmap(
    Final_matrix,
    name = "Membership",
    col = hm_cols,
    cluster_rows = F,
    cluster_columns = F,
    top_annotation = top_ann,
    column_split = col_split,
    column_title = NULL,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = TRUE,
    column_names_side = "top",
    row_title = NULL, 
    show_column_dend = F,
    row_split = cluster_factor,       
    left_annotation = cluster_ann, 
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(1, 2, 3),
      labels = c("Stat1 only", "Irf1 only", "Both"),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"), 
      labels_gp = gpar(fontsize = legend_size)    
    ),
    border = "black",
    width = unit(column_width, "mm") * ncol(Final_matrix)
  )
  
  draw(ht_final)
  
  
  return(
    list(
      heatmap = ht_final,
      rows = Cluster_df$gene,
      cols = unique(res$col_order),
      matrix = Final_matrix,
      cluster_df = Cluster_df
    )
  )
}


# ORDERING SLICES ---------------------------------

HM_orderSlices = function(mat){
  

  set.seed(12)
  
  # 3) Draw heatmap with the split and block annotation
  ht_1 = Heatmap(
    mat,
    name = "Membership",
    cluster_rows = TRUE,
    cluster_columns = F,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = TRUE,
    show_column_dend = F,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(1, 2, 3),
      labels = c("Stat1 only", "Irf1 only", "Both"),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"),  
      labels_gp = gpar(fontsize = legend_size)    
    ),
    border = "black"
  )
  
  ht_drawn = draw(ht_1)
  
  gene_order = row.names(mat)[row_order(ht_drawn)]
  
  return(gene_order)
}

# BUILD MEMBERSHIP BUCKETS ---------------------------------


build_membership_buckets = function(membership_df,
                                    celltype_priority = NULL,
                                    min_n = 20,
                                    type_priority = c("Both","Stat1_only","Irf1_only"),
                                    unique_celltype_assignment = FALSE) {
  stopifnot(all(c("Organ","CellType_Core","gene","has_Stat1","has_Irf1") %in% names(membership_df)))
  
  # cat per row
  df = membership_df %>%
    dplyr::mutate(
      cat = dplyr::case_when(
        has_Stat1 & has_Irf1 ~ "Both",
        has_Stat1 & !has_Irf1 ~ "Stat1_only",
        !has_Stat1 & has_Irf1 ~ "Irf1_only",
        TRUE ~ "None"
      )
    ) %>%
    dplyr::filter(cat != "None") %>%
    dplyr::mutate(
      CellType_Core = as.character(CellType_Core),
      gene = as.character(gene),
      cat = factor(cat, levels = c("Stat1_only","Irf1_only","Both"))
    ) %>%
    select(CellType_Core, gene, cat)
  
  df = df %>% dplyr::mutate(ct_rank = match(CellType_Core, celltype_priority))
  cat_levels = c("Both", "Stat1_only", "Irf1_only")
  
  
  df_unique =
    df %>%
    dplyr::arrange(ct_rank, CellType_Core, gene, match(cat, cat_levels))  %>%
    dplyr::group_by(gene, cat) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()%>%
    dplyr::arrange(ct_rank, CellType_Core, cat, gene)
  
  
  # 4) Buckets
  df_bucketed = df_unique %>%
    dplyr::transmute(CellType_Core,
                     orig_Type = as.character(cat),
                     Type = orig_Type,
                     gene,
                     bucket = paste(CellType_Core, Type, sep = "|"))
  
  bucket_list = df_bucketed %>%
    dplyr::group_by(bucket) %>%
    dplyr::summarise(genes = list(sort(unique(gene))), .groups = "drop") %>%
    { setNames(.$genes, .$bucket) }
  
  # ensure all buckets exist
  all_ct = sort(unique(df$CellType_Core))
  all_cat = c("Stat1_only","Irf1_only","Both")
  all_buckets = as.vector(outer(all_ct, all_cat, paste, sep = "|"))
  missing_b = setdiff(all_buckets, names(bucket_list))
  if (length(missing_b) > 0) {
    for (b in missing_b) bucket_list[[b]] = character(0)
    bucket_list = bucket_list[all_buckets]
  }
  
  # Cluster
  col_order = celltype_priority
  cluster_df = map_dfr(seq_along(bucket_list), function(i){
    
    cluster_name = names(bucket_list)[i]
    
    genes = unlist(bucket_list[i], use.names = FALSE)
    
    if (length(genes) == 0) return(NULL)
    
    data.frame(cluster_start = i,
               cluster_name = cluster_name, 
               gene = genes, 
               stringsAsFactors = FALSE)
  }) %>%
    separate(cluster_name, into = c("CellType_Core","Type"), sep = "\\|") %>%
    arrange(match(CellType_Core, col_order), match(Type, all_cat)) %>%
    mutate(cluster = dplyr::consecutive_id(cluster_start)) %>%
    select(cluster, CellType_Core, Type, gene) %>%
    left_join(
      .,
      df %>%
        select(CellType_Core, gene, cat) %>%
        distinct() %>%
        dplyr::rename(orig_Type = cat),
      by = c("CellType_Core","gene")
    )
  
  rownames(cluster_df) = NULL
  
  # merge small buckets
  counts0 = cluster_df %>%
    dplyr::count(CellType_Core, Type, name = "n") %>%
    tidyr::complete(CellType_Core = all_ct, Type = all_cat, fill = list(n = 0))
  
  map_merge =
    counts0 %>%
    dplyr::group_by(CellType_Core) %>%
    dplyr::mutate(
      Type_priority_rank = match(Type, type_priority),
      ord = order(dplyr::desc(n), Type_priority_rank),
      best_Type = Type[ord][1],
      Type_merge_to = dplyr::if_else(n < min_n, best_Type, Type)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(CellType_Core, Type, Type_merge_to)
  
  cluster_df_merged =
    cluster_df %>%
    dplyr::left_join(map_merge, by = c("CellType_Core","Type")) %>%
    dplyr::mutate(Type_merged = dplyr::coalesce(Type_merge_to, Type)) %>%
    dplyr::select(CellType_Core, Type = Type_merged, orig_Type, gene)
  
  # Recompute clusters 
  cluster_df_merged =
    cluster_df_merged %>%
    dplyr::mutate(bucket = paste(CellType_Core, Type, sep = "|")) %>%
    dplyr::arrange(match(CellType_Core, col_order), match(Type, all_cat), gene) %>%
    dplyr::group_by(CellType_Core, Type) %>%
    dplyr::mutate(cluster = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(match(CellType_Core, col_order), match(Type, all_cat)) %>%
    dplyr::mutate(cluster = dplyr::consecutive_id(paste(CellType_Core, Type))) %>%
    dplyr::select(cluster, CellType_Core, Type, orig_Type, gene)
  
  counts_final =
    cluster_df_merged %>%
    dplyr::count(cluster, CellType_Core, Type, name = "n") %>%
    dplyr::arrange(cluster)
  
  list(
    cluster_df = cluster_df_merged,        
    col_order = unique(col_order),
    counts = counts_final,
    params = list(min_n = min_n,
                  type_priority = type_priority,
                  unique_celltype_assignment = unique_celltype_assignment)
  )
}



# MEMBERSHIP HEATMAP BUCKETS ---------------------------------

Membership_HM_clusterlfc = function(matrix, col_order, cluster_df, core_df, reg_df){
  
  num_mat = matrix
  
  # Heatmap colors
  hm_cols =  c("0" = "#FFFFFF", "1" = "#B86A6A", "2" = "#5A7DA6", "3" = "#6EA27A") # none, Stat1, Irf1, Both
  
  # top annotation
  col_colors = core_palette[names(core_palette) %in% col_order]
  
  
  top_ann = HeatmapAnnotation(
    Organ = anno_block(
      gp = gpar(fill = col_colors, col = "black")
    ),
    height = unit(2, "mm")
  )
  
  
  Core_reference = Core_LFC_HM(
    core_df=core_df, 
    reg_df=reg_df, 
    row_order = cluster_df$gene, 
    col_order = col_order 
  )
  Core_reference_matrix = Core_reference$matrix
  Membership_reference_matrix = num_mat[cluster_df$gene,col_order]
  
  
  
  Cluster_df = map_dfr(unique(cluster_df$cluster), function(clust){
    
    df_slice = cluster_df %>% filter(cluster == clust)
    genes = df_slice$gene
    
    mat_slice_df = Membership_reference_matrix[genes,] %>%
      as.df %>%
      mutate(
        summary = do.call(paste, c(.[, col_order, drop = FALSE], sep = "|")),
        Type = consecutive_id(summary)
      )
    
    df_reorder = map_dfr(unique(mat_slice_df$Type), function(type){
      
      genes2 = mat_slice_df %>% filter(Type == type) %>% row.names
      
      if(length(genes2)<3){
        return( df_slice %>% filter(gene %in% genes2) )
      }
      
      mat_slice_lfc = Core_reference_matrix[genes2,]
      gene_order = HM_orderSlices(mat_slice_lfc)
      
      df_slice %>% filter(gene %in% genes2) %>% arrange(match(gene, gene_order))
      
    })
    
    
  })
  
  
  
  
  col_split = factor(col_order, levels = col_order)
  
  
  Final_matrix = num_mat[Cluster_df$gene,col_order]
  
  
  cluster_factor = factor(Cluster_df$cluster, levels = unique(Cluster_df$cluster))
  stopifnot(nrow(Final_matrix) == length(cluster_factor))
  
  
  # Cluster anno
  cluster_cols = generate_greys(
    max(as.numeric(Cluster_df$cluster))
  )
  
  cluster_ann = rowAnnotation(
    Cluster = anno_block(
      gp = gpar(fill = unname(cluster_cols), col = "black"),
      labels = levels(cluster_factor),
      labels_gp = gpar(col = "black", fontsize = 12)
    ),
    width = unit(5, "mm")
  )
  
  # Heatmap
  
  ht_final = Heatmap(
    Final_matrix,
    name = "Membership",
    col = hm_cols,
    cluster_rows = F,
    cluster_columns = F,
    top_annotation = top_ann,
    column_split = col_split,
    column_title = NULL,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = TRUE,
    column_names_side = "top",
    row_title = NULL, 
    show_column_dend = F,
    row_split = cluster_factor,       
    left_annotation = cluster_ann, 
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(1, 2, 3),
      labels = c("Stat1 only", "Irf1 only", "Both"),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"),   # legend title size
      labels_gp = gpar(fontsize = legend_size)     # legend element label size
    ),
    border = "black",
    width = unit(column_width, "mm") * ncol(Final_matrix)
  )
  
  draw(ht_final)
  
  
  return(
    list(
      heatmap = ht_final,
      rows = Cluster_df$gene,
      cols = col_order,
      matrix = Final_matrix,
      cluster_df = Cluster_df
    )
  )
}



# MEMBERSHIP_HM_ANNOTATED -------------------------------------------------------


Membership_HM_anno = function(Cluster_df, matrix, col_order, pathway_df, type_pathway = c("Wiki", "GO_BP"), pathway_count=10, pathway_w=NULL){
  
  Final_matrix = matrix[Cluster_df$gene,col_order]
  
    # Cluster annotations - left
  Cluster_df = Cluster_df %>% 
    group_by(cluster) %>%
    mutate(total = n()) %>%
    ungroup() %>%
    mutate(
      cluster_number = cluster,
      cluster = paste0(cluster, " (", total, ")")
    )
  
  cluster_factor = factor(Cluster_df$cluster, levels = unique(Cluster_df$cluster))
  stopifnot(nrow(Final_matrix %>% as.matrix) == length(cluster_factor))
  
  cluster_cols = generate_greys(
    max(seq_along(unique((Cluster_df$cluster_number))))
  )
  
  cluster_ann = rowAnnotation(
    Cluster2 = anno_block(
      gp = gpar(fill = unname(cluster_cols), col = "black"),
    ),
    width = unit(2, "mm")
  )
  
  
  # Pathway annotations
  
  pathway_palette_Go = c(
    "black",
    "#4d4d4d"
  )
  
  pathway_palette_wiki = c(
    "#172163",
    "#626999"
  )
  
  
  
  pathway_df = pathway_df %>%
    filter(source %in% type_pathway) %>%
    group_by(source, cluster) %>%
    arrange(Adjusted.P.value) %>%
    dplyr::slice(1:as.numeric(pathway_count)) %>%
    ungroup() %>%
    arrange(cluster, source, Adjusted.P.value)
  
  
  if(pathway_count >8){
    pathway_font = 6
  } else{pathway_font = 8}
  
  
  pathway_text_list = map(unique(pathway_df$cluster), function(clust){
    
    pathways = 
      rbind(
        data.frame(
          text = paste0("Cluster: ", clust),
          fontsize = as.numeric(pathway_font) +1,
          col = "#147d68"
        ),
        pathway_df %>%
          filter(cluster == clust) %>%
          mutate(
            whichColor = if_else(row_number() %% 2L == 1L, 1L, 2L),
            col = case_when(
              source == "Wiki" ~ pathway_palette_wiki[whichColor],
              source == "GO_BP" ~ pathway_palette_Go[whichColor],
              TRUE ~ "white"
            ),
            fontsize = pathway_font,
            text = Term
          ) %>%
          select(text, fontsize, col)
      )
    
    return(pathways)
    
  }) %>%
    set_names(levels(cluster_factor))
  
  
  
  pathway_anno = rowAnnotation(
    pathways = anno_textbox(
      pathway_text_list,
      which      = "row",
      align_to   = cluster_factor,
      side       = "left",
      background_gp = grid::gpar(fill = "#fafafa", col = "black"),
    )
  )
  
  if(!is.null(pathway_w)){
    pathway_anno = rowAnnotation(
      pathways = anno_textbox(
        pathway_text_list,
        which      = "row",
        align_to   = cluster_factor,
        side       = "left",
        background_gp = grid::gpar(fill = "#fafafa", col = "black"),
      ),
      width = max_text_width(pathway_text_list) + unit(pathway_w, "mm")
    )
  }
  
  cluster_ann_numbers = rowAnnotation(
    Cluster = anno_block(
      gp = gpar(fill = unname(cluster_cols), col = NA),
      labels = unique(Cluster_df$cluster),
      labels_gp = gpar(col = "black", fontsize = 8),
      labels_rot = 0 ,
      labels_just = "left"
    ),
    width = max_text_width(levels(cluster_factor), gp = gpar(fontsize = 11)) + unit(2, "mm")
  )
  
  lab_levels = levels(cluster_factor)
  lab_gp     = gpar(col = "black", fontsize = 8)
  
  cluster_ann_numbers = rowAnnotation(
    Cluster = anno_block(
      gp          = gpar(fill = unname(cluster_cols), col = NA),
      labels      = lab_levels,         
      labels_gp   = lab_gp,
      labels_rot  = 0,
      labels_just = "left",
      labels_offset = unit(2, "mm") 
    ),
    width = max_text_width(lab_levels, gp = lab_gp) + unit(4, "mm")
  )
  
  
  # Core annotation - top
  
  
  if(is.null(colnames(Final_matrix))){
    Final_matrix = as.matrix(Final_matrix)
    colnames(Final_matrix) = col_order
  }
  core_iterations = colnames(Final_matrix)
  
  core_factor = factor(core_iterations, levels = unique(core_iterations))
  stopifnot(ncol(Final_matrix) == length(core_factor))
  
  # distinct fills for each core cell type
  core_colors = data.frame(core_palette) %>%
    rownames_to_column("Core_CellTypes") %>%
    filter(Core_CellTypes %in% levels(core_factor)) %>%
    arrange(match(Core_CellTypes, unique(core_iterations))) %>%
    deframe()
  
  
  
  core_ann = HeatmapAnnotation(
    core_block = anno_block(
      gp = gpar(fill = core_colors, col = "black")
    ),
    which = "column",
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    height = unit(2, "mm")
  )
  
  
  
  ht_final =  Heatmap(
    Final_matrix,
    name = "Membership",
    col = c("0" = "white", "1" = "#c25555",
            "2" = "#5A7DA6", 
            "3" = "#6EA27A"),
    cluster_rows = F,
    cluster_columns = F,
    column_split = core_factor,
    top_annotation = core_ann,
    show_row_names = F,
    show_row_dend = F,
    show_column_names = TRUE,
    column_names_side = "top",
    row_title = NULL,
    column_title = NULL,
    show_column_dend = F,
    row_split = cluster_factor,
    left_annotation = c(pathway_anno, cluster_ann, cluster_ann_numbers),
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    show_heatmap_legend = T,
    heatmap_legend_param = list(
      at = c(1, 2, 3),
      labels = c("Stat1 only", "Irf1 only", "Both"),
      title_gp  = gpar(fontsize = title_size, fontface = "bold"),   # legend title size
      labels_gp = gpar(fontsize = legend_size)     # legend element label size
    ),
    border = "black",
    width = unit(column_width, "mm") * ncol(Final_matrix)
  )
  
  
  return(
    list(
      heatmap = ht_final,
      rows = Cluster_df$gene,
      cols = col_order,
      matrix = Final_matrix,
      cluster_df = Cluster_df
    )
  )
}









