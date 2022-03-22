
# library(synExtra)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(here)
library(gt)
library(qs)
library(fst)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data", .cache = TRUE)

paste_ <- function(...) {
  paste(..., sep = "_")
}

theme_set(theme_light())

wd <- here("fig1")
dir.create(wd, showWarning = FALSE)

raw_counts <- syn("syn25292308") %>%
  qread()

meta <- syn("syn25292310") %>%
  qread()

perturbation_meta <- syn("syn21547097") %>%
  read_csv()

MIN_N_READS_DETECTABLE <- 5
MIN_N_READS_VALID <- 50000

sequencing_stats <- raw_counts %>%
  rowwise() %>%
  transmute(
    dataset, plate, date,
    stats = counts %>%
      group_by(sample_id) %>%
      summarize(
        total_counts = sum(count),
        n_detectable = sum(count >= MIN_N_READS_DETECTABLE),
        fraction_detectable = n_detectable / n(),
        passing = total_counts >= MIN_N_READS_VALID
      ) %>%
      list()
  ) %>%
  ungroup() %>%
  unnest(stats)

dataset_names <- tribble(
  ~dataset, ~dataset_name, ~date,
  "jsm_merck", "Myofibroblast screen", "2017_09",
  "ld_dub", "DUB inhibitors", "2018_06",
  "sr_repurposing", "AD drug repurposing 1", "2018_02",
  "sr_repurposing", "AD drug repurposing 2", "2018_11",
  "lincs_cdk4_6_7", "CDK4/6 inhibitors 1", "2017_09",
  "lincs_cdk4_6_7", "CDK4/6 inhibitors 2", "2019_08",
  "fp_transdiff", "2015_10_feodor_price_transdifferentiation_screen", "2015_10",
  "okl", "Optimal Kinase Library", "2021_08"
) %>%
  mutate(across(c(dataset_name), fct_inorder))

seq_depth_violin <- sequencing_stats %>%
  filter(!dataset %in% c("fp_transdiff")) %>%
  inner_join(dataset_names, by = c("dataset", "date")) %>%
  mutate(across(dataset_name, fct_rev)) %>%
  ggplot((aes(x = total_counts, y = dataset_name))) +
    geom_quasirandom(varwidth = FALSE, width = 0.45, method = "quasirandom", groupOnX = FALSE) +
    scale_x_log10() +
    scale_y_discrete(position = "left") +
    geom_vline(xintercept = MIN_N_READS_VALID, linetype = "dashed") +
    labs(y = "", x = "Sequencing depth")
    # theme(
    #   # axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
    #   # axis.text.y = element_text(face = "bold")
    # )

ggsave(
  file.path(wd, "sequencing_depth_violin.pdf"),
  seq_depth_violin,
  width = 5.5, height = 4.5
)

detectable_genes_violin <- sequencing_stats %>%
  filter(!dataset %in% c("fp_transdiff"), passing) %>%
  inner_join(dataset_names, by = c("dataset", "date")) %>%
  mutate(across(dataset_name, fct_rev)) %>%
  ggplot((aes(x = n_detectable, y = dataset_name, color = passing))) +
    geom_quasirandom(varwidth = FALSE, width = 0.45, method = "quasirandom", groupOnX = FALSE) +
    scale_x_log10() +
    scale_y_discrete(position = "left") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red"), guide = FALSE) +
    # geom_vline(xintercept = MIN_N_READS_VALID, linetype = "dashed") +
    labs(y = "", x = "N genes with >5 reads per library")


ggsave(
  file.path(wd, "detectable_genes_violin.pdf"),
  detectable_genes_violin, width = 5.5, height = 4.5
)

cor_depth_detectable_plot <- sequencing_stats %>%
  filter(!dataset %in% c("fp_transdiff"), passing) %>%
  inner_join(dataset_names, by = c("dataset", "date")) %>%
  ggplot((aes(x = total_counts, y = n_detectable, color = passing))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red"), guide = FALSE) +
    facet_wrap(vars(dataset_name))
    scale_x_discrete(position = "top") +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "red"), guide = FALSE) +
    # geom_hline(yintercept = MIN_N_READS_VALID, linetype = "dashed") +
    labs(x = "", y = "N genes with >5 reads") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
      axis.text.y = element_text(face = "bold")
    )

dataset_names <- tribble(
  ~dataset, ~dataset_name, ~date,
  "jsm_merck", "1) Myofibroblast screen", "2017_09",
  "ld_dub", "2) DUB inhibitors [xx]", "2018_06",
  "lincs_cdk4_6_7", "3) CDK4/6 inhibitors 1 [xx]", "2017_09",
  "lincs_cdk4_6_7", "4) CDK4/6 inhibitors 2 [xx]", "2019_08",
  "sr_repurposing", "5) AD drug repurposing 1 [xx]", "2018_02",
  "sr_repurposing", "6) AD drug repurposing 2 [xx]", "2018_11",
  "fp_transdiff", "2015_10_feodor_price_transdifferentiation_screen", "2015_10",
  "okl", "7) Optimal Kinase Library [xx]", "2021_08"
) %>%
  mutate(across(c(dataset_name), fct_inorder))

meta_table <- meta %>%
  filter(!dataset %in% c("fp_transdiff")) %>%
  inner_join(dataset_names, by = c("dataset", "date")) %>%
  unnest(meta) %>%
  filter(
    # Don't want this cell line from Laura
    !cells %in% c("MDAMB231"),
    # In CDK dataset, only use the four CDK4/6 inhibs, not CDK7
    if_else(
      dataset == "lincs_cdk4_6_7",
      lspci_id %in% c(
        89588, 91464, 78621, 94539
      ) | drug == "control",
      TRUE
    )
  ) %>%
  group_by(dataset, dataset_name, date) %>%
  summarize(
    `n samples` = n(),
    `n drugs` = length(unique(na.omit(drug_id))),
    `n drugs in CMap` = length(intersect(lspci_id, perturbation_meta$lspci_id)),
    `n cell lines` = length(unique(cells)),
    .groups = "drop"
  ) %>%
  inner_join(
    sequencing_stats %>%
      select(dataset, date, total_counts, n_detectable) %>%
      group_by(dataset, date) %>%
      summarize(
        `median sequencing depth` = median(total_counts) %>%
          as.integer(),
        `median genes detected` = median(n_detectable) %>%
          as.integer()
      ),
    by = c("dataset", "date")
  ) %>%
  select(-dataset, -date) %>%
  arrange(dataset_name) %>%
  rename(` ` = dataset_name)

meta_table_gt <- meta_table %>%
  gt() %>%
  tab_options(column_labels.font.weight = "bold") %>%
  cols_align("right") %>%
  cols_align("left", columns = vars(` `)) %>%
  fmt_number(
    columns = tidyselect::eval_select(rlang::expr(where(is.numeric)), meta_table),
    use_seps = TRUE,
    sep_mark = ",",
    decimals = 0
  )

gtsave(
  meta_table_gt,
  file.path(wd, "fig1b.html")
)


