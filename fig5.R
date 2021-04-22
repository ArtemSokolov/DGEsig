library(synExtra)
library(tidyverse)
library(cmapR)
library(here)
library(cowplot)
library(broom)
library(ggrepel)
library(seriation)

synapser::synLogin()
syn <- synExtra::synDownloader("data")

wd <- here("fig5")
dir.create(wd, showWarnings = FALSE)

theme_set(theme_bw())

compound_name_map <- syn("syn22035396.3") %>%
  read_rds() %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1)

cmap_gene_meta <- syn("syn21547102") %>%
  read_csv()

clue_res_dge <- syn("syn21907139") %>%
  read_rds()

clue_res_l1000 <- syn("syn21907143") %>%
  read_rds()

clue_res_combined <- syn("syn21907166.4") %>%
  read_rds()
# clue_res_combined_2 <- clue_res_combined

pertubation_meta <- syn("syn21547097.6") %>%
  read_csv()

signature_meta <- syn("syn21547101.6") %>%
  read_csv()

dge_meta <- syn("syn22000707.8") %>%
  read_rds()

deseq_res <- syn("syn22017733.1") %>%
  read_rds()

deseq_meta <- syn("syn21558154.4") %>%
  read_rds()

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14),
        axis.ticks = element_blank())
}


R <- clue_res_combined %>%
  filter( result_type == "pert", score_level == "summary" ) %>%
  pluck( "data", 1 ) %>%
  # Fix compound names
  mutate(
    name_query = if_else(
      str_count(name_query, "[0-9]") > 1 | str_detect(name_query, fixed("-")),
      str_to_upper(name_query),
      name_query
    )
  )

## Identify relevant gene sets
relevant_gene_sets <- deseq_res %>%
  mutate(
    gene_set = condition_conc %>%
      str_replace_all("\\s", "_") %>%
      str_replace_all("[^\\w]", "")
  ) %>%
  inner_join(
    dge_meta %>%
      unnest(meta) %>%
      mutate(
        query_group = case_when(
          lspci_id %in% c(90768, 91047, 93750, 77650) ~ "controls",
          !lspci_id %in% signature_meta$lspci_id ~ paste("unknowns", dataset, sep = " "),
          TRUE ~ NA_character_
        )
      ) %>%
      filter(
        !is.na(query_group),
        dataset %in% c("sr_repurposing", "ld_dub", "lincs_cdk4_6_7"),
        is.na(stim)
      ) %>%
      group_by(query_group, dataset, cells, drug_id) %>%
      arrange(desc(time)) %>%
      slice(1) %>%
      select(query_group, cells, drug_id, time, stim, dataset),
    by = c("cells", "drug_id", "time", "stim")
  )

## List cell lines used for profiling compounds
cell_lines_used <- relevant_gene_sets %>%
  group_by(lspci_id) %>%
  summarize(
    cells = case_when(
      all(
        c(
          "BT549", "HCC1806", "Hs578T", "MCF7", "PDX1258", "PDXHCI002", "T47D"
        ) %in% cells
      ) ~ if ("rencell" %in% cells) "breast cancer & RenCells"
          else "8 breast cancer cell lines",
      cells == "rencell" ~ "differentiated RenCells",
      TRUE ~ cells
    ), .groups = "drop"
  ) %>%
  distinct() %>%
  inner_join(
    distinct(R, name_query, lspci_id_query),
    by = c("lspci_id" = "lspci_id_query")
  )

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- R %>%
  drop_na(tau) %>%
  # filter(pert_type == "trt_cp") %>%
  inner_join(
    relevant_gene_sets %>%
      distinct(query_group, gene_set, drug_id, dataset),
    by = "gene_set"
  ) %>%
  # dplyr::rename(pert_iname = id) %>%
  group_by( query_group, pert_type, pert_iname, name_query, source, z_score_cutoff ) %>%
  summarize_at(
    "tau", ~quantile(.x, c(0.67, 0.33), names = FALSE) %>%
      {.[order(abs(.))[2]]}
    # "tau", ~.x[ which.max(abs(.x)) ]
  ) %>%
  ungroup() %>%
  mutate_at( "source", toupper )

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
      group_by(., name_query) %>%
        arrange(desc(abs(tau)), .by_group = TRUE) %>%
        slice(1:10) %>%
        ungroup() %>%
        pull(pert_iname)
    }
  )


## Perform hierarchical clustering on drugT profiles (columns in the final plot)
## Use the DGE slice because it is cleaner and less saturated
DM <- R3 %>%
  filter( source == "DGE" ) %>%
  select(
    name_query,
    pert_iname,
    tau
  ) %>%
  spread( name_query, tau ) %>%
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
lvl2 <- R2 %>%
  filter( source == "DGE" ) %>%
  distinct(name_query, query_group) %>%
  group_by(query_group) %>%
  group_map(
    function(.x, ...) {
      queries <- unique(.x[["name_query"]])
      if (length(queries) == 1)
        return(queries)
      DM[, queries] %>%
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
      }
  }

cell_lines_used_plot_data <- cell_lines_used %>%
  bind_rows(
    crossing(
      name_query = setdiff(lvl2, unique(R3[["name_query"]])),
      y = "a",
      cells = NA_character_
    )
  ) %>%
  crossing(split_group = c(1, 2)) %>%
  mutate(
    name_query = factor(name_query, levels = lvl2),
    y = "a"
  )

cell_lines_used_plot <- cell_lines_used_plot_data %>%
  # Arrange just so that cells with border are drawn last for clean rendering
  arrange(!str_detect(name_query, "^( )+$")) %>%
  ggplot(aes(x = name_query, y = y, fill = cells)) +
    geom_tile(aes(color = str_detect(name_query, "^( )+$"))) +
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
      guide = FALSE
    ) +
    # Remove facet labels
    theme(
      strip.background = element_blank(), strip.text.x = element_blank(),
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.text.y = element_blank(), axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

## Fix the order via factor levels
R4 <- R3 %>%
  bind_rows(
    crossing(
      name_query = setdiff(lvl2, unique(R3[["name_query"]])),
      pert_iname = unique(.[["pert_iname"]]),
      tau = NA_real_
    )
  ) %>%
  mutate(
    # Horizontal facet split unite the second and third dendrogram cut
    split_group = split_vector[pert_iname] %>%
      as.character(),
    # recode("2" = "3"),
    name_query = factor(name_query, lvl2),
    pert_iname = factor(pert_iname, lvl)
  )

## Plotting a heatmap of clue hits
fplot <- function(X, cell_lines = NULL) {
  ggplot(
    # Arrange just so that cells with border are drawn last for clean rendering
    X %>% arrange(!str_detect(name_query, "^( )+$")),
    aes(x=name_query,
        y=pert_iname,
        fill=tau)
  )+
    theme_minimal() + theme_bold() +
    geom_tile(
      aes(color = str_detect(name_query, "^( )+$"))
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
      panel.grid = element_blank()
    )
}

connectivity_plot <- fplot(R4)

combined_plot <- plot_grid(
  cell_lines_used_plot,
  connectivity_plot,
  ncol = 1,
  align = "v", axis = "lr",
  rel_heights = c(1, 25)
)

ggsave(
  file.path(wd, "fig5.pdf"),
  combined_plot,
  width = 11.5, height = 10
)

