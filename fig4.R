library(synExtra)
library(tidyverse)
library(cmapR)
library(here)
library(egg)
library(broom)
library(ggrepel)
library(fst)
library(qs)
library(data.table)

synapser::synLogin()
syn <- synExtra::synDownloader("data", .cache = TRUE)

wd <- here("fig4")
dir.create(wd, showWarnings = FALSE)

theme_set(theme_bw())

pertubation_meta <- syn("syn21547097") %>%
  fread()

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  bind_rows(
    anti_join(pertubation_meta, ., by = "lspci_id") %>%
      select(name = pert_iname, lspci_id) %>%
      drop_na(name)
  ) %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()

cmap_gene_meta <- syn("syn21547102") %>%
  fread()

diff_exp_by_conc <- syn("syn25303743") %>%
  fread()

diff_exp_by_conc_agg <- diff_exp_by_conc %>%
  filter(
    concentration_method == "per_concentration", replicate_method == "replicates_aggregated",
    stim == "" | is.na(stim)
  )
  # # Don't aggregate time points
  # group_by(
  #   drug_id, lspci_id, cells, drug_conc, ensembl_gene_id
  # ) %>%
  # summarize(
  #   # across(
  #   #   c(log2FoldChange, log2FoldChange_MLE),
  #   #   ~quantile(.x, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
  #   #     {.[order(abs(.))[2]]},
  #   # ),
  #   across(
  #     c(log2FoldChange, log2FoldChange_MLE),
  #     mean
  #   ),
  #   .groups = "drop"
  # )

signature_meta <- syn("syn21547101") %>%
  fread()

cmap_signatures_mcf7 <- syn("syn27254561") %>%
  qread() %>%
  filter(cell_aggregate_method == "per_cell_line", replicate_method == "replicates_aggregated") %>%
  chuck("data", 1) %>%
  inner_join(
    distinct(cmap_gene_meta, pr_gene_id, entrez_id, ensembl_gene_id, symbol) %>%
      drop_na()
  )

ensembl_hugo_map <- diff_exp_by_conc %>%
  pull(ensembl_gene_id) %>%
  unique() %>%
  genebabel::query_hgnc("ensembl_gene_id")

etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
        axis.text.y = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14),
        title = etxt(12))
}

calculate_dose_response_l1000 <- function(signatures, z_threshold = 1.645) {
  signatures %>%
    group_by(entrez_id) %>%
    filter(any(abs(zscore) > z_threshold)) %>%
    ungroup() %>%
    mutate(change_norm = zscore) %>%
    rename(change = zscore, drug_conc = pert_dose) %>%
    group_nest(entrez_id)
}

calculate_dose_response_dge <- function(data, p_threshold = 0.05) {
  # browser()
  selected_signatures <- data %>%
    # select(condition, result) %>%
    # unnest(result) %>%
    # mutate(log2FoldChange = log2FoldChange_MLE) %>%
    # group_by(condition) %>%
    mutate(
      change_norm = scale(log2FoldChange),
      change = log2FoldChange
    ) %>%
    # ungroup() %>%
    group_by(ensembl_gene_id) %>%
    filter(any(na.omit(padj) <= p_threshold)) %>%
    group_nest()
  selected_signatures
}

calculate_dose_response <- function(drug, cell, time, p_threshold = 0.05, method = "spearman", round_digits = 0) {
  l1000 <- cmap_signatures_mcf7 %>%
    filter(lspci_id == drug, cell_id == cell, if (!is.null(time)) pert_time == time else TRUE) %>%
    select_if(negate(is.list)) %>%
    calculate_dose_response_l1000(z_threshold = qnorm(1 - p_threshold)) %>%
    inner_join(cmap_gene_meta %>% filter(pr_is_lm == 1) %>% select(entrez_id, symbol, ensembl_gene_id), by = "entrez_id")
  dge <- diff_exp_by_conc_agg %>%
    filter(lspci_id == !!drug, cells == cell, if (!is.null(time)) replace_na(time == !!time, TRUE) else TRUE) %>%
    calculate_dose_response_dge(p_threshold = p_threshold) %>%
    inner_join(select(ensembl_hugo_map, ensembl_gene_id, symbol), by = "ensembl_gene_id")
  # browser()
  list(
    l1000 = l1000,
    dge = dge
  ) %>%
    bind_rows(.id = "method") %>%
    mutate(
      test = map(
        data,
        ~suppressWarnings(cor.test(
          if (round_digits > 0) round(.x[["change"]], digits = round_digits) else .x[["change"]],
          .x[["drug_conc"]], method = !!method
        ))
      )
    ) %>%
    group_by(symbol) %>%
    filter(length(unique(method)) == 2, n() == 2) %>%
    ungroup() %>%
    select(symbol, method, test, data)
}

dose_response_pairwise <- function(df) {
  df %>%
    select(symbol, method, test) %>%
    spread(method, test) %>%
    mutate(
      significance = select(., dge, l1000) %>%
        mutate_all(map_lgl, ~.x[["p.value"]] < 0.05) %>%
        {
          case_when(
            .[["dge"]] & .[["l1000"]] ~ "both",
            .[["dge"]] ~ "dge_only",
            .[["l1000"]] ~ "l1000_only",
            TRUE ~ "neither"
          )
        }
    ) %>%
    mutate_at(vars(dge, l1000), map_dbl, "estimate")
}

dose_response_pairwise_fisher <- function(df) {
  tab <- df %>%
    select(dge, l1000) %>%
    gather("method", "correlation") %>%
    mutate_at(vars(correlation), ~if_else(abs(.x) > 0.25, "correlated", "uncorrelated")) %>%
    table()
  test <- fisher.test(tab, conf.int = TRUE,conf.level = .95)
  list(tab, test)
}

dose_response_curve <- function(df) {
  scale_factor <- df %>%
    group_by(method) %>%
    summarize(max_val = max(abs(change))) %>%
    ungroup() %>%
    spread(method, max_val) %>%
    {.[["l1000"]]/.[["dge"]]}
  df %>%
    mutate(
      change_scaled = if_else(method == "dge", change*scale_factor, change)
    ) %>%
    # ggplot(aes(drug_conc, change_norm, color = method)) +
    ggplot(aes(drug_conc, change_scaled, color = method)) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_point(alpha = 0.5) +
    scale_y_continuous(
      sec.axis = sec_axis(~./scale_factor, name = "DGE")
    ) +
    geom_smooth(aes(fill = method, group = method), method = "lm", alpha = 0.1) +
    scale_x_log10() +
    labs(x = "Dose (uM)", y = "L1000")
}

dose_response_cor_and_curves <- function(df, seed = 42, highlighted_genes = NULL) {
  # browser()
  set.seed(seed)
  dose_cor_plot_data <- dose_response_pairwise(df) %>%
    mutate(
      pos = paste(
        case_when(
          dge < -0.25 ~ "left",
          dge < 0.25 ~ "middle",
          dge < 1 ~ "right",
          TRUE ~ NA_character_
        ),
        case_when(
          l1000 < -0.25 ~ "bottom",
          l1000 < 0.25 ~ "middle",
          l1000 < 1 ~ "top",
          TRUE ~ NA_character_
        ),
        sep = "_"
      )
    )
  # Select random gene from each position unless
  # we explicitly selected one
  highlighted_genes_all <- dose_cor_plot_data  %>%
    # Only interested in these sectors atm
    filter(
      pos %in% c("left_middle", "left_bottom", "middle_middle", "right_top", "right_middle", "middle_top", "middle_bottom")
    ) %>%
    group_by(pos) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    # Remove positions where a gene has been manually selected
    anti_join(
      dose_cor_plot_data %>%
        filter(symbol %in% highlighted_genes),
      by = "pos"
    ) %>%
    pull(symbol) %>%
    c(highlighted_genes)
  dose_cor_plot_data <- dose_cor_plot_data %>%
    mutate(
      selected = if_else(symbol %in% highlighted_genes_all, pos, "no")
    )
  dose_fisher <- dose_response_pairwise_fisher(dose_cor_plot_data)
  cor_plot <- dose_cor_plot_data %>%
    mutate(
      color_selected = if_else(selected == "no", "black", "red") %>%
        factor(levels = c("black", "red"))
    ) %>%
    arrange(color_selected) %>%
    ggplot(aes(dge, l1000, color = color_selected)) +
    geom_hline(yintercept = c(-0.25, 0.25), color = "grey30") +
    geom_vline(xintercept = c(-0.25, 0.25), color = "grey30") +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "yellow", alpha = 0.2)),
      xmin = -1, xmax = -0.25, ymin = -1, ymax = -0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "yellow", alpha = 0.2)),
      xmin = 0.25, xmax = 1, ymin = 0.25, ymax = 1
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "magenta", alpha = 0.2)),
      xmin = -0.25, xmax = 0.25, ymin = 0.25, ymax = 1
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "magenta", alpha = 0.2)),
      xmin = -0.25, xmax = 0.25, ymin = -1, ymax = -0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "cyan", alpha = 0.2)),
      xmin = -1, xmax = -0.25, ymin = -0.25, ymax = 0.25
    ) +
    annotation_custom(
      grid::rectGrob(gp = grid::gpar(fill = "cyan", alpha = 0.2)),
      xmin = 0.25, xmax = 1, ymin = -0.25, ymax = 0.25
    ) +
    # coord_cartesian(expand = FALSE) +
    geom_point(aes(size = selected, text = symbol), alpha = 0.8) +
    scale_size_manual(
      # values = c("no" = 2, "left_middle" = 3, "left_bottom" = 3, "right_middle" = 3, "right_top" = 3),
      values = c(set_names(rep_len(3, length.out = length(unique(dose_cor_plot_data$pos))), unique(dose_cor_plot_data$pos)), "no" = 2),
      guide = FALSE
    ) +
    geom_text_repel(
      aes(label = symbol),
      data = ~.x %>%
        mutate(symbol = if_else(selected == "no", "", symbol)),
      box.padding = 0.5,
      point.padding = 0.5,
      max.iter = 8000
    ) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    scale_color_identity() +
    guides(color = FALSE) +
    labs(x = "Dose-response correlation DGE", y = "Dose-respose correlation L1000")
  dose_plots <- df %>%
    inner_join(
      dose_cor_plot_data %>%
        filter(selected != "no"),
      by = "symbol"
    ) %>%
    select(symbol, pos, method, data) %>%
    mutate(data = map(data, select, -symbol)) %>%
    unnest(data) %>%
    group_nest(symbol, pos) %>%
    mutate(
      plot = map(data, dose_response_curve) %>%
        map2(symbol, ~.x + labs(title = .y) + theme(aspect.ratio = 0.5))
    ) %>%
    {set_names(.[["plot"]], .[["pos"]])}
  # browser()
  res <- arrangeGrob(
    grobs = list(
      cor_plot +
        labs(title = "Palbociclib dose-response correlation"),
      dose_plots[["left_middle"]] +
        theme(legend.position = "none"),
      dose_plots[["left_bottom"]] +
        theme(legend.position = "none"),
      dose_plots[["right_middle"]],
      dose_plots[["right_top"]],
      dose_plots[["middle_top"]] +
        theme(legend.position = "none"),
      dose_plots[["middle_middle"]] +
        theme(legend.position = "none"),
      dose_plots[["middle_bottom"]]
    ) %>%
      map(~.x + theme_bold()) %>%
      c(
        list(
          grid::nullGrob(),
          grid::nullGrob()
        )
      ),
    layout_matrix = rbind(
      c(2, 9, 1, 10, 5),
      c(3, 9, 1, 10, 4),
      c(6, 9, 7, 10, 8)
    ),
    widths = unit(c(4, 0.1, 6, 0.1, 5), "in"),
    heights = unit(c(3, 3, 3), "in")
  )
  attr(res, "center_gg") <- cor_plot
  res
}

palbo_dose_cor <- calculate_dose_response(
  filter(compound_names, name == "PALBOCICLIB")$lspci_id,
  "MCF7", time = 24, method = "pearson", round_digits = 0,
  p_threshold = 0.05
)

picked_genes <- c(
  "CDC25B", "HACD3", "DUSP4", "CTSD"
)
# For selecting genes
# attr(palbo_dose_cor_plot, "center_gg") %>% ggplotly()

palbo_dose_cor_plot <- dose_response_cor_and_curves(
  palbo_dose_cor, seed = 42, highlighted_genes = picked_genes
)

grid::grid.draw(palbo_dose_cor_plot)

ggsave(
  file.path(wd, "palbo_dose_cor_plot.pdf"),
  palbo_dose_cor_plot, width = 15.5, height = 6
)

palbo_dose_cor_plot_data <- dose_response_pairwise(palbo_dose_cor)

palbo_dose_fisher <- dose_response_pairwise_fisher(palbo_dose_cor_plot_data)

palbo_dose_cor_plot <- palbo_dose_cor_plot_data %>%
  ggplot(aes(dge, l1000)) +
  geom_hline(yintercept = c(-0.25, 0.25), color = "grey30") +
  geom_vline(xintercept = c(-0.25, 0.25), color = "grey30") +
  geom_point(alpha = 0.8) +
  coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(title = "Palbociclib dose-response correlation per gene")

ggsave(
  file.path(wd, "palbo_dose_cor_plot_significance.pdf"),
  palbo_dose_cor_plot,
  width = 6, height = 5
)


palbo_dose_srsf3 <- palbo_dose_cor %>%
  select(symbol, method, data) %>%
  filter(symbol == "SRSF3") %>%
  unnest(data) %>%
  dose_response_curve()

ggsave(
  file.path(wd, "palbo_dose_srsf3.pdf"),
  palbo_dose_srsf3,
  width = 6, height = 4
)

# Calculate dose-response correlation for all drugs where we have MCF-7 data
# and run a Fisher's test to determine if more genes show dose-response
# relationship with DGE or with L1000

dose_response_all <- diff_exp_by_conc_agg %>%
  distinct(lspci_id, cells, time) %>%
  inner_join(
    select(cmap_signatures_mcf7, lspci_id, cells = cell_id, time = pert_time)
  ) %>%
  distinct() %>%
  mutate(
    data = pmap(
      list(lspci_id, cells, time),
      possibly(calculate_dose_response, NULL), method = "pearson", round_digits = 0,
      p_threshold = 0.05
    )
  ) %>%
  filter(!map_lgl(data, is.null)) %>%
  mutate(
    pairwise_data = map(
      data, dose_response_pairwise
    ),
    fisher_res = map(
      pairwise_data, dose_response_pairwise_fisher
    ),
    fisher_tibble = map(fisher_res, 2) %>%
      map(broom::tidy)
  ) %>%
  unnest(fisher_tibble)

dose_response_all %>% select(-data, -pairwise_data, -fisher_res) %>% View()

p <- ggplot(
  dose_response_all, aes(-log10(p.value))
) +
  geom_histogram()

dose_response_all_plotting <- dose_response_all %>%
  mutate(
    signed_p = -log10(p.value) * if_else(estimate > 0, 1, -1),
    ordered_lspci_id = paste(lspci_id, cells, time) %>%
      factor(levels = .[order(signed_p)])
  )

p <- ggplot(
  dose_response_all_plotting,
  aes(
    ordered_lspci_id, signed_p, fill = ordered_lspci_id == paste(filter(compound_names, name == "PALBOCICLIB")$lspci_id, "MCF7", "24")
  )
) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "grey50"), guide = "none") +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() +
  theme_bold() +
  theme(axis.text.x = element_blank()) +
  labs(x = "Rank ordered drugs", y = "-log10(p)")
  # coord_cartesian(expand = FALSE)

ggsave(
  file.path(wd, "all_drugs_correlation_pvalues.pdf"),
  p,
  width = 5, height = 3
)

