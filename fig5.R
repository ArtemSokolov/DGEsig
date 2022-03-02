library(synExtra)
library(tidyverse)
library(cmapR)
library(here)
library(cowplot)
library(broom)
library(ggrepel)
library(seriation)
library(qs)

synapser::synLogin()
syn <- synExtra::synDownloader("data")

wd <- here("fig5")
dir.create(wd, showWarnings = FALSE)

theme_set(theme_bw())

cmap_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = cols(
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  bind_rows(
    anti_join(cmap_meta, ., by = "lspci_id") %>%
      select(name = pert_iname, lspci_id) %>%
      drop_na(name)
  ) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()


cmap_gene_sets <- syn("syn25314203") %>%
  qread() %>%
  as_tibble()

dge_gene_sets <- syn("syn25303778") %>%
  qread() %>%
  as_tibble()

gene_set_meta <- bind_rows(
  l1000 = cmap_gene_sets %>%
    select(
      cell_aggregate_method, replicate_method,
      lspci_id, drug_conc = pert_dose, time = pert_time, cells = cell_id,
      cutoff, gene_set_id
    ),
  dge = dge_gene_sets %>%
    select(
      concentration_method, replicate_method,
      lspci_id, drug_conc, cells, stim, stim_conc,
      time, gene_set_id
    ),
  .id = "source"
) %>%
  mutate(cells = coalesce(cells, "summary"))

x <- inner_join(
  dge_gene_sets,
  cmap_gene_sets %>%
    filter(cell_aggregate_method == "per_cell_line"),
  by = c("lspci_id", "cells" = "cell_id")
) %>%
  distinct(lspci_id, cells)

cmap_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = cols(
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

sig_meta <- syn("syn21547101") %>%
  read_csv()

dge_meta <- syn("syn25292310") %>%
  qread()

x <- inner_join(
  dge_meta %>%
    unnest(meta),
  sig_meta,
  by = c("lspci_id", "cells" = "cell_id")
) %>%
  distinct(lspci_id, cells) %>%
  drop_na() %>%
  left_join(
    compound_names
  )


pertubation_meta <- syn("syn21547097") %>%
  read_csv()

signature_meta <- syn("syn21547101") %>%
  read_csv()

# dge_meta <- syn("syn25292310") %>%
#   qread()

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14),
        axis.ticks = element_blank())
}

## Load all results
R <- syn("syn26468923") %>%
  # file.path(pathData, "clue_results_combined.rds") %>%
  qread() %>%
  filter( result_type == "pert", score_level == "summary" ) %>%
  pluck( "data", 1 )

## Identify relevant gene sets
relevant_gene_sets <- dge_gene_sets %>%
  drop_na(lspci_id) %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    lspci_id %in% {dge_meta %>%
      unnest(meta) %>%
      filter(
          !lspci_id %in% cmap_gene_sets$lspci_id,
          dataset %in% c("sr_repurposing", "ld_dub", "lincs_cdk4_6_7", "okl")
      ) %>%
      pull(lspci_id) %>%
        c(94539, 91464, 78621, 89588)},
    is.na(stim)
  ) %>%
  mutate(
    query_group = case_when(
      lspci_id %in% c(94539, 91464, 78621, 89588) ~ "controls",
      TRUE ~ "unknowns"
    )
  )
  # filter(
  #   concentration_method == "concentration_aggregated",
  #   replicate_method == "replicates_aggregated",
  #   !lspci_id %in% cmap_gene_sets$lspci_id,
  #   dataset %in% c("sr_repurposing", "ld_dub", "lincs_cdk4_6_7"),
  #   is.na(stim)
  # ) %>%


## List cell lines used for profiling compounds
cell_lines_used <- relevant_gene_sets %>%
  distinct(lspci_id, cells) %>%
  drop_na() %>%
  group_by(lspci_id) %>%
  summarize(
    cells = if(
      all(
        c(
          "BT549", "HCC1806", "Hs578T", "MCF7", "PDX1258", "PDXHCI002", "T47D"
        ) %in% cells
      )
    ) if ("rencell" %in% cells) "breast cancer & RenCells"
          else "8 breast cancer cell lines"
    else if (all(cells == "rencell")) "differentiated RenCells"
    else paste(cells, collapse = " "),
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = "name")
  )

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- R %>%
  drop_na(tau) %>%
  # filter(pert_type == "trt_cp") %>%
  inner_join(
    relevant_gene_sets %>%
      distinct(query_group, gene_set_id, lspci_id),
    by = c("gene_set" = "gene_set_id")
  ) %>%
  # dplyr::rename(pert_iname = id) %>%
  group_by( query_group, pert_type, pert_iname, lspci_id ) %>%
  summarize_at(
    "tau", ~quantile(.x, c(0.67, 0.33), names = FALSE) %>%
      {.[order(abs(.))[2]]}
    # "tau", ~.x[ which.max(abs(.x)) ]
  ) %>%
  ungroup()

# FInd most similar drugs per query

R_only_genetic <- R2 %>%
  filter(pert_type %in% c("trt_oe", "trt_sh.cgs")) %>%
  mutate(
    pert_iname = paste(pert_iname, pert_type, sep = "_") %>%
      str_replace(fixed("_trt_oe"), " OE") %>%
      str_replace(fixed("_trt_sh.cgs"), " KD")
  )

# Filter to only include pertubatios that are in the top 10 of any query
R3 <- R_only_genetic %>%
  dplyr::filter(
    pert_iname %in% {
      group_by(., lspci_id) %>%
        arrange(desc(abs(tau)), .by_group = TRUE) %>%
        slice(1:10) %>%
        ungroup() %>%
        pull(pert_iname)
    }
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = "name")
  )


## Perform hierarchical clustering on drugT profiles (columns in the final plot)
## Use the DGE slice because it is cleaner and less saturated
DM <- R3 %>%
  select(
    name_q,
    pert_iname,
    tau
  ) %>%
  spread( name_q, tau ) %>%
  as.data.frame %>%
  column_to_rownames("pert_iname")

comp_order <- function(clust) {
  dm <- clust %>%
    dist()
  dm %>%
    hclust() %>%
    reorder(dm) %>% {
      .[["labels"]][.[["order"]]]
    }
}

DM_clust <- DM %>%
  dist() %>%
  hclust() %>%
  reorder(dist(DM))

lvl <- DM_clust %>% {
  .[["labels"]][.[["order"]]]
}

# Divide heatmap into two pieces that we can plot side-by-side
split_vector <- set_names(
  rep(c(1, 2), each = ceiling(length(lvl) / 2), length.out = length(lvl)),
  lvl
)

# Complicated way to fix factor levels on x-axis so it includes spaces
# between groups of query datasets
lvl2 <- R3 %>%
  distinct(name_q, query_group) %>%
  group_by(query_group) %>%
  group_map(
    function(.x, ...) {
      queries <- unique(.x[["name_q"]])
      if (length(queries) == 1)
        return(queries)
      DM[, as.character(queries)] %>%
        t() %>%
        comp_order()
    }
  ) %>%
  # Add blank spaces between datasets
  {
    imap(
      .,
      ~c(.x, strrep(" ", .y))
    ) %>%
      {
        exec(c, !!!.)
      } %>%
      head(n = -1)
  }

cell_lines_used_plot_data <- cell_lines_used %>%
  bind_rows(
    crossing(
      name_q = setdiff(lvl2, unique(R3[["name_q"]])),
      y = "a",
      cells = NA_character_
    )
  ) %>%
  crossing(split_group = c(1, 2)) %>%
  mutate(
    name_q = factor(name_q, levels = lvl2),
    y = "a"
  )

cell_lines_used_plot <- cell_lines_used_plot_data %>%
  # Arrange just so that cells with border are drawn last for clean rendering
  arrange(!str_detect(name_q, "^( )+$")) %>%
  ggplot(aes(x = name_q, y = y, fill = cells)) +
    geom_tile(aes(color = str_detect(name_q, "^( )+$"))) +
    facet_wrap(~split_group, nrow = 1) +
    theme_minimal() + theme_bold() +
    scale_fill_brewer(
      palette = "Set2",
      breaks = na.omit(cell_lines_used_plot_data$cells),
      na.value = "white"
    ) +
    # Remove cell borders in empty space between left and right side
    scale_color_manual(
      values = c("TRUE" = NA_character_, "FALSE" = "black"),
      guide = "none"
    ) +
    # Remove facet labels
    theme(
      strip.background = element_blank(), strip.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(), axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

## Fix the order via factor levels
R4 <- R3 %>%
  bind_rows(
    crossing(
      name_q = setdiff(lvl2, unique(R3[["name_q"]])),
      pert_iname = unique(.[["pert_iname"]]),
      tau = NA_real_
    )
  ) %>%
  mutate(
    # Horizontal facet split unite the second and third dendrogram cut
    split_group = split_vector[pert_iname] %>%
      as.character(),
    # recode("2" = "3"),
    name_q = factor(name_q, lvl2),
    pert_iname = factor(pert_iname, lvl)
  )

## Plotting a heatmap of clue hits
fplot <- function(X, cell_lines = NULL) {
  ggplot(
    # Arrange just so that cells with border are drawn last for clean rendering
    X %>% arrange(!str_detect(name_q, "^( )+$")),
    aes(x=name_q,
        y=pert_iname,
        fill=tau)
  )+
    theme_minimal() + theme_bold() +
    geom_tile(
      aes(color = str_detect(name_q, "^( )+$"))
    ) +
    # Remove cell borders in empty space between left and right side
    scale_color_manual(
      values = c("TRUE" = NA_character_, "FALSE" = "black"),
      guide = FALSE
    ) +
    scale_fill_gradientn( colors=pal, limits=c(-100,100), na.value = "white" ) +
    labs(x = "Drug query", y = "Clue target class", fill = "Tau" ) +
    facet_wrap(~split_group, scales = "free_y") +
    # Remove facet labels
    theme(
      strip.background = element_blank(), strip.text.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.text.y = element_blank(), axis.title.y = element_blank()
    )
}

connectivity_plot <- fplot(R4)

combined_plot <- plot_grid(
  connectivity_plot,
  cell_lines_used_plot,
  ncol = 1,
  align = "v", axis = "lr",
  rel_heights = c(25, 10)
)

ggsave(
  file.path(wd, "fig5.pdf"),
  combined_plot,
  width = 11.5, height = 10
)

