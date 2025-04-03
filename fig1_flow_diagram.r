library(tidyverse)
library(synExtra)
library(here)
library(qs)
library(powerjoin)

pathData <- "~/data"

wd <- here("fig1")
dir.create(wd, showWarnings = FALSE)

synapser::synLogin()
syn <- synDownloader(pathData, .cache = TRUE)

raw_counts <- syn("syn25292308") %>%
  qread()

perturbation_meta <- syn("syn21547097") %>%
  read_csv()

## Load all results
R_all <- syn("syn26468923") %>%
  # file.path(pathData, "clue_results_combined.rds") %>%
  qread()

R <- R_all %>%
  # filter( result_type == "pert", score_level == "cell" ) %>%
  # pluck( "data", 1 ) %>%
  # filter(cell_id == "MCF7")
  filter( result_type == "pert", score_level == "summary" ) %>%
  pluck( "data", 1 ) %>%
  as_tibble()

M <- syn("syn25292310") %>%
  # syn("syn22000707") %>%
  qread() %>%
  unnest(meta)

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

meta_filtered <- M %>%
  filter(!dataset %in% c("fp_transdiff")) %>%
  inner_join(dataset_names, by = c("dataset", "date")) %>%
  filter(
    sample_id %in% {
      raw_counts$counts %>%
        map("sample_id") %>%
        map(unique) %>%
        {rlang::exec("c", !!!.)}
    },
    # Don't want this cell line from Laura
    !cells %in% c("MDAMB231"),
    # Exclude non-drug rencell treatments
    !drug_id %in% c("dsrnalipo", "lps", "nakeddsrna", "lipocontrol"),
    # In CDK dataset, only use the four CDK4/6 inhibs, not CDK7
    if_else(
      dataset == "lincs_cdk4_6_7",
      lspci_id %in% c(
        89588, 91464, 78621, 94539
      ) | drug == "control",
      TRUE
    ),
    is.na(stim) | stim == "control"
  )

meta_filtered %>%
  distinct()

cmap_gene_sets <- syn("syn25314203.4") %>%
  qread() %>%
  filter(
    cell_aggregate_method == "cells_aggregated",
    replicate_method == "replicates_aggregated"
  )

dge_gene_sets <- syn("syn25303778") %>%
  qread() %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    is.na(stim) | stim == "control"
  ) %>%
  as_tibble() %>%
  semi_join(
    meta_filtered,
    by = c("drug_id", "cells", "time")
  )

meta_conditions <- meta_filtered %>%
  filter(drug_norm != "control") %>%
  distinct(
    cells, time, drug_id, lspci_id
  )

meta_conditions %>%
  anti_join(dge_gene_sets)

dge_gene_sets %>%
  anti_join(meta_conditions)

nrow(meta_conditions)
nrow(dge_gene_sets)

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
  unnest(stats) %>%
  semi_join(
    meta_filtered
  )

nrow(sequencing_stats)

sequencing_stats %>%
  count(passing)
# # A tibble: 2 × 2
#   passing     n
#   <lgl>   <int>
# 1 FALSE     175
# 2 TRUE     2760

meta_conditions %>%
  count(cells)
# 9 cell lines

meta_conditions %>%
  count(drug_id)
# 133 drugs

meta_conditions %>%
  nrow()
# 230 conditions

nrow(dge_gene_sets)
# 227 conditions

dge_gene_sets %>%
  count(drug_id)
# 130 drugs

dge_gene_sets %>%
  count(cells)
# 9 cell lines

R_oi <- R %>%
  distinct(gene_set) %>%
  semi_join(
    dge_gene_sets %>%
      select(gene_set_id) %>%
      distinct(),
    by = c("gene_set" = "gene_set_id")
  )



dge_gene_sets %>%
  anti_join(
    R,
    by = c("gene_set_id" = "gene_set")
  ) %>%
  select(-gene_sets, -ends_with("method")) %>%
  print(n = Inf)

gs_w_R <- dge_gene_sets %>%
  semi_join(
    R,
    by = c("gene_set_id" = "gene_set")
  ) %>%
  select(-gene_sets, -ends_with("method"))

gs_w_R %>%
  count(cells)
# 9 cell lines

gs_w_R %>%
  count(drug_id)
# 102 drugs

nrow(gs_w_R)
# 172 conditions

R %>%
  distinct(gene_set) %>%
  anti_join(
    dge_gene_sets,
    by = c("gene_set" = "gene_set_id")
  )

R %>%
  distinct(gene_set) %>%
  semi_join(
    dge_gene_sets,
    by = c("gene_set" = "gene_set_id")
  )



meta_cell_lines <- tribble(
  ~dataset, ~date, ~cells,
  "jsm_merck", "2017_09", "Myofibroblasts",
  "ld_dub", "2018_06", "MCF7",
  "lincs_cdk4_6_7", "2017_09", "7 breast cancer lines",
  "lincs_cdk4_6_7", "2019_08", "3 breast cancer lines",
  "okl", "2021_08", "MCF7",
  "sr_repurposing", "2018_02", "ReNcell VM",
  "sr_repurposing", "2018_11", "ReNcell VM"
)

library(gt)

meta_table <- meta_filtered %>%
  filter(drug_norm != "control") %>%
  group_by(dataset, dataset_name, date) %>%
  summarize(
    `n samples` = n(),
    `n drugs` = length(unique(na.omit(drug_id))),
    `n drugs in CMap` = length(intersect(lspci_id, perturbation_meta$lspci_id)),
    # `n cell lines` = length(unique(cells)),
    .groups = "drop"
  ) %>%
  power_inner_join(
    meta_cell_lines %>%
      rename(`Cell lines` = cells),
    by = c("dataset", "date"),
    check = check_specs(
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn",
      duplicate_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  power_inner_join(
    sequencing_stats %>%
      select(dataset, date, total_counts, n_detectable) %>%
      group_by(dataset, date) %>%
      summarize(
        `median sequencing depth` = signif(median(total_counts) / 1e5, 2),
        `median genes detected` = signif(median(n_detectable) / 1e3, 2)
      ),
    by = c("dataset", "date"),
    check = check_specs(
      unmatched_keys_left = "warn",
      unmatched_keys_right = "warn",
      duplicate_keys_left = "warn",
      duplicate_keys_right = "warn"
    )
  ) %>%
  select(-dataset, -date) %>%
  arrange(dataset_name) %>%
  rename(` ` = dataset_name)

meta_table_gt <- meta_table %>%
  gt() %>%
  opt_table_font("Helvetica") %>%
  tab_options(column_labels.font.weight = "bold") %>%
  cols_align("right") %>%
  cols_align("left", columns = vars(` `)) %>%
  cols_label(
    `median sequencing depth` ~ "Median sequencing depth x{{10^5}}",
    `median genes detected` ~ "Median genes detected x{{10^3}}"
  )
  # fmt_number(
  #   columns = tidyselect::eval_select(rlang::expr(where(is.numeric)), meta_table),
  #   use_seps = TRUE,
  #   sep_mark = ",",
  #   decimals = 0
  # )

gtsave(
  meta_table_gt,
  file.path(wd, "fig1b.html")
)
