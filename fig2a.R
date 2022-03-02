
# library(synExtra)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(here)
library(gt)
library(ComplexUpset)
library(ggvenn)
library(qs)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data", .cache = TRUE)

paste_ <- function(...) {
  paste(..., sep = "_")
}

theme_set(theme_light())

wd <- here("fig2")
dir.create(wd, showWarning = FALSE)

meta <- syn("syn25292310") %>%
  qread() %>%
  filter(dataset != "fp_transdiff") %>%
  unnest(meta)

perturbation_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = cols(
      lspci_id = col_integer(),
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )



cmap_gene_sets <- syn("syn25314203") %>%
  qread()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  bind_rows(
    anti_join(perturbation_meta, ., by = "lspci_id") %>%
      select(name = pert_iname, lspci_id) %>%
      drop_na(name)
  ) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()


dge_all_compounds <- meta %>%
  filter(dataset != "fp_transdiff") %>%
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
  distinct(dataset, lspci_id) %>%
  drop_na() %>%
  left_join(
    compound_names
  )

good_dge_sets <- dge_gene_sets %>%
  filter(concentration_method == "concentration_aggregated", replicate_method == "replicates_aggregated") %>%
  unnest(gene_set_table) %>%
  filter(padj < 0.1) %>%
  group_by(lspci_id) %>%
  filter(
    if(sum(direction == "up") >= 10 && sum(direction == "down") >= 10) TRUE else FALSE
  ) %>%
  ungroup()

clue_res_all <- syn("syn26468923") %>%
  qread()

clue_res_summary <- clue_res_all %>%
  filter(
    result_type == "pert",
    score_level == "summary"
  ) %>%
  chuck("data", 1) %>%
  left_join(
    perturbation_meta %>%
      distinct(pert_id, lspci_id),
    by = "pert_id"
  )


upset_data <- meta %>%
  distinct(lspci_id) %>%
  drop_na() %>%
  mutate(
    dge = TRUE
  ) %>%
  left_join(
    perturbation_meta %>%
      filter(dataset != "LINCS_2020") %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        cmap = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    perturbation_meta %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        cmap_new = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    good_dge_sets %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(
        gene_set = TRUE
      ),
    by = "lspci_id"
  ) %>%
  left_join(
    clue_res_summary %>%
      distinct(lspci_id) %>%
      drop_na() %>%
      mutate(cmap_returns = TRUE),
    by = "lspci_id"
  ) %>%
  mutate(
    across(everything(), replace_na, replace = FALSE)
  )

upset_plot <- upset(
  upset_data,
  intersect = c("dge", "cmap", "cmap_new", "gene_set", "cmap_returns"),
  set_sizes = upset_set_size(
    geom = geom_bar()
  ) +
    geom_text(
      aes(label=stat(count)),
      stat='count',
      color = "white",
      hjust = -0.5
    )
)

cowplot::ggsave2(
  file.path(wd, "fig2a.pdf"),
  upset_plot, width = 6, height = 4
)
