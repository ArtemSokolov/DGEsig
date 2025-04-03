library(tidyverse)
library(powerjoin)
# library(DESeq2)
library(synapser)
library(here)
library(qs)
library(data.table)

synLogin()

syn <- synExtra::synDownloader("~/data", .cache = TRUE)
# syn <- synExtra::synDownloader("data", .cache = TRUE)

dge_gene_sets <- syn("syn25303778") %>%
  qread()

cmap_gene_sets <- syn("syn27768305") %>%
  qread()

cmap_perturbation_meta <- syn("syn21547097") %>%
  read_csv(
    col_types = list(
      cmap_name = col_character(),
      compound_aliases = col_character(),
      target = col_character(),
      moa = col_character()
    )
  )

cmap_signature_meta <- syn("syn21547101") %>%
  read_csv()

cmap_gene_meta <- syn("syn21547102") %>%
  fread()

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  bind_rows(
    anti_join(cmap_perturbation_meta, ., by = "lspci_id") %>%
      select(name = pert_iname, lspci_id) %>%
      drop_na(name)
  ) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()


gene_set_comparison <- dge_gene_sets %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "per_replicate",
    is.na(stim),
    !is.na(lspci_id)
  ) %>%
  transmute(
    technique = "dge",
    lspci_id, cells,
    replicate = as.factor(paste(dataset, time, replicate, sep = "_")),
    time,
    gene_set = gene_set_table %>%
      map(~transmute(.x, direction, score=-log10(padj), entrezgene_id))
  ) %>%
  bind_rows(
    cmap_gene_sets %>%
      filter(
        replicate_method == "per_replicate",
        cell_aggregate_method == "per_cell_line",
        cutoff == 0.7
      ) %>%
      chuck("data", 1) %>%
      semi_join(
        dge_gene_sets,
        by = c("lspci_id", "cell_id" = "cells")
      ) %>%
      transmute(
        technique = "l1000", lspci_id, cells = cell_id, replicate,
        gene_set = map(
          gene_set_table,
          \(x) power_inner_join(
            x,
            select(cmap_gene_meta, pr_gene_id, entrezgene_id = entrez_id),
            by = "pr_gene_id",
            check = check_specs(
              unmatched_keys_left = "warn",
              duplicate_keys_right = "warn"
            )
          ) %>%
            transmute(direction, score = abs(zscore), entrezgene_id)
        )
      ) %>%
      power_inner_join(
        cmap_signature_meta %>%
          distinct(replicate = sig_id, time = pert_time),
        by = "replicate",
        check = check_specs(
          unmatched_keys_left = "warn",
          duplicate_keys_right = "warn"
        )
      )
  ) %>%
  arrange(technique, desc(time)) %>%
  group_by(technique, lspci_id, cells, time) %>%
  mutate(
    name = paste0(str_to_upper(technique), " ", time, "h Rep ", rev(seq_len(n()))) %>%
      fct_inorder()
  ) %>%
  ungroup()

gene_set_comparison_venn_data <- gene_set_comparison %>%
  mutate(across(cells, str_to_upper)) %>%
  # mutate(gene_set = map(gene_set, group_nest, direction, .key = "gene_set")) %>%
  # unnest(gene_set) %>%
  # nest(gene_set = c(direction, score, entrezgene_id)) %>%
  group_by(lspci_id, cells) %>%
  filter(n() > 1) %>%
  summarize(
    # venn = cur_data() %>%
    venn_data = cur_data() %>%
      mutate(across(gene_set, map, arrange, desc(score))) %>%
      list(),
    .groups = "drop"
  )



make_upset <- function(df, ...) {
  df %>%
    transmute(
      name,
      gene_set = map(
        gene_set, head, n = 300
      )
    ) %>%
    unnest(gene_set) %>%
    mutate(dummy = TRUE) %>%
    distinct(entrezgene_id, name, dummy) %>%
    pivot_wider(
      id_cols = entrezgene_id,
      names_from = name,
      values_from = dummy,
      values_fill = FALSE
    ) %>%
    select(-entrezgene_id)
    # UpSetR::upset(
    #   nsets = ncol(.),
    #   keep.order = FALSE,
    #   order.by = "freq",
    #   group.by = "degree",
    #   ...
    # )
}

set.seed(3)
gene_set_comparison_upset_data <- gene_set_comparison_venn_data %>%
  left_join(
    cmap_perturbation_meta %>%
      distinct(lspci_id, pert_iname) %>%
      drop_na(),
    by = "lspci_id"
  ) %>%
  filter(
    # make sure that there are at least two replicates for dge and l1000
    map_lgl(venn_data, ~{tbl <- table(.x[["technique"]]); length(tbl) >= 2 && all(tbl >= 2)})
  ) %>%
  mutate(
    # Selecting 5 random L1000 samples at 24h
    venn_data = map(
      venn_data,
      ~.x %>%
        group_by(technique) %>%
        filter(
          if (technique[1] == "dge")
            TRUE
          else {
            time == 24
          }
        ) %>%
        filter(
          if (technique[1] == "dge")
            TRUE
          else {
            seq_len(n()) %in% sample(seq_len(n()), 5)
          }
        )
    ),
    upset_data = map(venn_data, make_upset)
  )

all_combn <- function(x) {
  map(
    seq_along(x),
    ~combn(x, .x, simplify = FALSE)
  ) %>%
    reduce(c)
}

gene_set_comparison_upset <- gene_set_comparison_upset_data %>%
  # head(n = 3) %>%
  mutate(
    venn = map(
      upset_data,
      ~{
        # browser()
        ComplexUpset::upset(
          .x,
          intersect = colnames(.x),
          n_intersections = 30,
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size()
            # 'Intersection ratio' = ComplexUpset::intersection_ratio()
          ),
          # queries = c(
          #   map(
          #     filter(.x, str_starts(name, "DGE")) %>%
          #       pull(name) %>%
          #       all_combn(),
          #     ~{
          #       # browser()
          #       ComplexUpset::upset_query(
          #         intersect = as.character(.x),
          #         only_components = "intersections_matrix",
          #         color = "cornflowerblue",
          #         fill = "cornflowerblue"
          #       )
          #     }
          #   ),
          #   map(
          #     filter(.x, str_starts(name, "L1000")) %>%
          #       pull(name) %>%
          #       all_combn(),
          #     ~ComplexUpset::upset_query(
          #       intersect = as.character(.x),
          #       only_components = "intersections_matrix",
          #       color = "orangered",
          #       fill = "orangered"
          #     )
          #   )
          # ),
          set_sizes = FALSE,
          sort_sets = FALSE,
          min_degree = 2,
          height_ratio = 2
        )
      }
    )
  )

pwalk(
  gene_set_comparison_upset,
  function(pert_iname, venn, cells, ...) {
    withr::with_pdf(
      paste0("fig3b_upset_de_genes_", pert_iname, "_", cells, ".pdf"),
      width = 12, height = 4,
      print(venn)
    )
  }
)

single_set_pie_charts <- gene_set_comparison_upset %>%
  mutate(
    pie = map(
      upset_data,
      function(x) {
        single_set <- x %>%
          as.matrix() %>%
          rowSums() %>%
          magrittr::is_greater_than(1)
        df <- data.frame(
          group = c("Single gene set", "Multiple gene sets"),
          value = c(sum(!single_set), sum(single_set))
        ) %>%
          arrange(desc(group)) %>%
          mutate(prop = value / sum(value)) %>%
          mutate(ypos = cumsum(value)- 0.5*value )
        ggplot(df, aes(x="", y=value, fill=group)) +
          geom_bar(stat="identity", width=1, color="white") +
          coord_polar("y", start=0) +
          labs(fill = "") +
          geom_text(aes(y = ypos, label = scales::percent(prop, 1)), color = "white", size=6) +
          theme_void() # remove background, grid, numeric labels
      }
    )
  )

pwalk(
  single_set_pie_charts,
  function(pert_iname, pie, cells, ...) {
    ggsave(
      paste0("fig3b_single_set_pie_", pert_iname, "_", cells, ".pdf"),
      pie, width = 3, height = 2.5
    )
  }
)

