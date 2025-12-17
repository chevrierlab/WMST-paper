# Setup -----------------------------------------------------------

source("E:/Coding/Scripts/R_Sourcing/R_Master_Source.R")


# Load data --------------------------------------------------------

CellType_Deconvolution_df = read.csv2("Z:/nchevrier/projects/clevenger1/Projects/ArraySeq/General_Analysis/Revision_Analysis/Benchmark_deconvolution/Annotations/deconvolution_results.csv")


# Generate agreement plot ------------------------------------------

organs = c("LI", "SP", "KI", "LU")
Decon_methods = c("Card", "RCTD", "Cell2Location")

sample_df = CellType_Deconvolution_df %>%
  filter(Decon_Method != "CellKB", Organ %in% organs) %>%
  group_by(Decon_Method, Organ, Sample) %>%
  summarize(correct = mean(Match == "Yes") * 100, .groups = "drop") %>% 
  arrange(match(Organ, organs), match(Decon_Method, Decon_methods), Sample) %>%
  mutate(Organ = factor(Organ, levels = organs), Decon_Method = factor(Decon_Method, c("Card", "RCTD", "Cell2Location"))) %>%
  as.df

means_df = sample_df %>%
  group_by(Decon_Method, Organ) %>%
  summarize(correct = mean(correct), .groups = "drop") %>% 
  arrange(match(Organ, organs), match(Decon_Method, Decon_methods)) %>%
  mutate(Organ = factor(Organ, levels = organs), Decon_Method = factor(Decon_Method, c("Card", "RCTD", "Cell2Location"))) %>%
  as.df

summary_plot = ggplot(means_df, aes(x = Organ, y = correct, fill = Decon_Method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(
    data = sample_df,
    aes(x = Organ, y = correct, group = Decon_Method),
    color = "black",
    position = position_jitterdodge(
      dodge.width  = 0.8,
      jitter.width = 0.08,
      jitter.height = 0
    ),
    size = 1.5,
    shape = 16,
    stroke = 0.6,
    alpha = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Card" = "#3e9b7c",
      "RCTD" = "#7e7db2",
      "Cell2Location" = "#d86416"
    )
  ) +
  labs(
    x = "",
    y = "Spot-level agreement with CellKb annotations",
    fill = "Deconvolution Method"
  ) +
  expand_limits(y = c(0, 100)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "right"
  ) 

