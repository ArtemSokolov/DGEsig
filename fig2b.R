library( tidyverse )
library( seriation )   # For optimal leaf reordering
library(synapser)
library(here)
library(qs)
library(magrittr)
library(data.table)

pathData <- "~/data"

wd <- here("fig2")
dir.create(wd, showWarnings = FALSE)

synapser::synLogin()
syn <- synExtra::synDownloader(pathData, .cache = TRUE)

rdbu_colormap <- RColorBrewer::brewer.pal(5, "RdBu")

cmap_nonlinear_colormap <- list(
    colors = c(
        rdbu_colormap[1:3],
        rdbu_colormap[3:5]
    ) %>%
        rev(),
    values = scales::rescale(c(-100, -90, -80, 80, 90, 100), from = c(-100, 100))
)

## Load all results
R_all <- syn("syn26468923") %>%
    # file.path(pathData, "clue_results_combined.rds") %>%
    qread()

R <- R_all %>%
  # filter( result_type == "pert", score_level == "cell" ) %>%
  # pluck( "data", 1 ) %>%
  # filter(cell_id == "MCF7")
    filter( result_type == "pert", score_level == "summary" ) %>%
    pluck( "data", 1 )

condition_conc_vars <- c("cells", "drug_id", "lspci_id", "stim", "stim_conc", "time")

M <- syn("syn25292310") %>%
    # syn("syn22000707") %>%
    qread() %>%
    unnest(meta)

cmap_gene_sets <- syn("syn25314203.4") %>%
  qread()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

R_cmap_mcf7 <- syn("syn30410740") %>%
  fread() %>%
  filter(cell_id == "MCF7")

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

M_cell_info <- bind_rows(
  select(cmap_gene_sets, gene_set_id, lspci_id, cells = cell_id),
  select(dge_gene_sets, gene_set_id, lspci_id, cells)
) %>%
  distinct()

# pertubation_meta <- syn("syn21547097") %>%
#     read_csv()

# pertubation_meta <- syn("syn21547097.6") %>%
#     read_csv()

## Identify the common set of drugs between DGE-query, L1000-query and targets
cmap_returned <- filter(cmap_meta, pert_id %in% R$pert_id)$lspci_id %>% unique()
dge_queried <- filter(dge_gene_sets, gene_set_id %in% R$gene_set)$lspci_id %>% unique()
qcom <- intersect(cmap_returned, dge_queried) %>%
  na.omit()

dmap <- cmap_meta %>%
  filter(pert_id %in% R$pert_id) %>%
  pull(lspci_id) %>%
  intersect(qcom)

diff <- setdiff(qcom, dmap)

## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>% {
    bind_rows(
      filter(., source == "l1000"),
      # Actually didn't perform cell aggregation for dge data... need to do here
      filter(., source == "dge") %>%
        group_by(source, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_t) %>%
        summarize(
          tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
            {.[order(abs(.))[2]]},
          # tau = mean(tau),
          cells_q = paste(sort(unique(cells_q)), collapse = "_"),
          .groups = "drop"
        ) %>%
        mutate(cell_aggregate_method = "cells_aggregated")
    )
  } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, cell_aggregate_method, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_q, cells_t, cutoff) %>%
  summarize(
      tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
        {.[order(abs(.))[2]]},
    # tau = mean(tau),
      .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  left_join(
    compound_names %>%
      rename(name_t = name),
    by = c("lspci_id_t" = "lspci_id")
  )

R2_by_cell <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "per_cell_line" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    cutoff == 0.7 | is.na(cutoff),
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R_all %>%
      filter( result_type == "pert", score_level == "cell" ) %>%
      pluck( "data", 1 ) %>%
      filter(cell_id == "MCF7") %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, cell_aggregate_method, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_q, cells_t, cutoff) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    # tau = mean(tau),
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  left_join(
    compound_names %>%
      rename(name_t = name),
    by = c("lspci_id_t" = "lspci_id")
  )


## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2_mcf7_cmap <- R_cmap_mcf7 %>%
  mutate(
    cells_q = "MCF7"
  ) %>%
  rename(cells_t = cell_id) %>%
  select(starts_with("lspci_id"), starts_with("cells"), tau) %>%
  drop_na(starts_with("lspci_id")) %>%
  group_by(across(-tau)) %>%
  summarize(tau = mean(tau), .groups = "drop") %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  left_join(
    compound_names %>%
      rename(name_t = name),
    by = c("lspci_id_t" = "lspci_id")
  )

cluster_mat <- function(df) {
  # Start off with data where replicates, cells, and concentrations aggregated and
  # using a z-threshold of 0.7
  lspci_id_common <- intersect(df$lspci_id_t, df$lspci_id_q)
  R3 <- df %>%
    filter(
      lspci_id_t %in% lspci_id_common,
      lspci_id_q %in% lspci_id_common
    ) %>%
    mutate(
      name_q = str_trunc(name_q, 16, ellipsis = "…"),
      name_t = str_trunc(name_t, 16, ellipsis = "…")
    )

  ## Perform hierarchical clustering on drugT profiles (columns in the final plot)
  ## Use the DGE slice because it is cleaner and less saturated
  R3_mat <- R3 %>% filter(source == "dge") %>% select( name_q, name_t, tau ) %>%
    # IKK16 slipped through. It was profiled with DGE but
    # it's signature for whatever reason didn't return results in CMap
    # (too few genes?). L1000 signature was succesfull
    filter(name_q != "IKK16", name_t != "IKK16") %>%
    pivot_wider( names_from = name_t, values_from = tau ) %>%
    column_to_rownames("name_q") %>% as.matrix()
  # We want to cluster both rows and cols but want them ordered
  # identically. Just adding distance matrices of matrix and it's transpose
  # DM <- R3_mat %>% t() %>% dist()
  DM <- (R3_mat %>% t() %>% dist()) +
    R3_mat %>% dist()
  lvl <- hclust(DM) %>% reorder(DM) %>%  dendextend::order.hclust() %>% labels(DM)[.]

  ## Fix the order via factor levels
  R4 <- R3 %>% mutate(name_q = factor(name_q, lvl),
                      name_t = factor(name_t, rev(lvl)))

  # Complete missing observations at z-scores that yielded insufficient
  # genes for Clue with NA
  R4_completed <- bind_rows(
    R4 %>%
      filter(source == "dge"),
    R4 %>%
      filter(source == "l1000") %>%
      complete(nesting(lspci_id_q, lspci_id_t, name_q, name_t, source), cutoff)
  ) %>%
    # IKK16 slipped through. It was profiled with DGE but
    # it's signature for whatever reason didn't return results in CMap
    # (too few genes?). L1000 signature was succesfull
    filter(name_q != "IKK16", name_t != "IKK16")
  R4_completed
}

R4_neuro <- R2_by_cell %>%
  filter(cells_q == "rencell") %>%
  cluster_mat()

p_neuro_query <- fplot(R4_neuro)
ggsave(
  "fig2_neuro_query.pdf", p_neuro_query, width = 8, height = 7
)

R4_mcf7 <- R2_by_cell %>%
  filter(cells_q == "MCF7") %>%
  cluster_mat()

p_mcf7_query <- fplot(R4_mcf7)
ggsave(
  "fig2_mcf7_query.pdf", p_mcf7_query, width = 7, height = 6
)

R4_mcf7_cmap <- R2_mcf7_cmap %>%
  # just a workaround it's actually l1000
  mutate(source = "dge", cutoff = 0.7) %>%
  cluster_mat()

p_mcf7_cmap <- fplot(R4_mcf7_cmap)
ggsave(
  "fig2_mcf7_cmap_query.pdf", p_mcf7_cmap, width = 14, height = 12
)

neuro_mcf7_overlap <- intersect(
  R2_by_cell %>%
    filter(cells_q == "MCF7") %>%
    pull(lspci_id_q),
  R2_by_cell %>%
    filter(cells_q == "rencell") %>%
    pull(lspci_id_q)
)
R4_neuro_both <- R2_by_cell %>%
  filter(cells_q == "rencell", lspci_id_q %in% neuro_mcf7_overlap) %>%
  cluster_mat()
p_neuro_query_both <- fplot(R4_neuro_both)
ggsave(
  "fig2_neuro_query_both.pdf", p_neuro_query_both, width = 6, height = 5
)
R4_mcf7_both <- R2_by_cell %>%
  filter(cells_q == "MCF7", lspci_id_q %in% neuro_mcf7_overlap) %>%
  cluster_mat()

p_mcf7_query_both <- fplot(R4_mcf7_both)
ggsave(
  "fig2_mcf7_query_both.pdf", p_mcf7_query_both, width = 6, height = 5
)

# Start off with data where replicates, cells, and concentrations aggregated and
# using a z-threshold of 0.7
R3 <- R2 %>%
  filter(
    cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    replicate_method == "replicates_aggregated" | is.na(replicate_method),
    concentration_method == "concentration_aggregated" | is.na(concentration_method),
    # cutoff == 0.7 | is.na(cutoff),
    lspci_id_t %in% qcom,
    lspci_id_q %in% qcom
  ) %>%
  mutate(
    name_q = str_trunc(name_q, 16, ellipsis = "…"),
    name_t = str_trunc(name_t, 16, ellipsis = "…")
  )

## Perform hierarchical clustering on drugT profiles (columns in the final plot)
## Use the DGE slice because it is cleaner and less saturated
R3_mat <- R3 %>% filter(source == "dge") %>% select( name_q, name_t, tau ) %>%
  # IKK16 slipped through. It was profiled with DGE but
  # it's signature for whatever reason didn't return results in CMap
  # (too few genes?). L1000 signature was succesfull
  filter(name_q != "IKK16", name_t != "IKK16") %>%
  pivot_wider( names_from = name_t, values_from = tau ) %>%
  column_to_rownames("name_q") %>% as.matrix()
# We want to cluster both rows and cols but want them ordered
# identically. Just adding distance matrices of matrix and it's transpose
# DM <- R3_mat %>% t() %>% dist()
DM <- (R3_mat %>% t() %>% dist()) +
  R3_mat %>% dist()
lvl <- hclust(DM) %>% reorder(DM) %>%  dendextend::order.hclust() %>% labels(DM)[.]

## Fix the order via factor levels
R4 <- R3 %>% mutate(name_q = factor(name_q, lvl),
                    name_t = factor(name_t, rev(lvl)))

# Complete missing observations at z-scores that yielded insufficient
# genes for Clue with NA
R4_completed <- bind_rows(
    R4 %>%
        filter(source == "dge"),
    R4 %>%
        filter(source == "l1000") %>%
        complete(nesting(lspci_id_q, lspci_id_t, name_q, name_t, source), cutoff)
) %>%
  # IKK16 slipped through. It was profiled with DGE but
  # it's signature for whatever reason didn't return results in CMap
  # (too few genes?). L1000 signature was succesfull
  filter(name_q != "IKK16", name_t != "IKK16")

## Plotting elements
pal <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
    theme(axis.text.x = etxt(12, angle=90, hjust=1, vjust=0.5),
          axis.text.y = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14),
          axis.ticks = element_blank())
}

## Plotting a heatmap of clue hits
fplot <- function(X) {
    ggplot( X, aes(x=name_t, y=name_q, fill=tau) ) +
        theme_minimal() + theme_bold() +
        geom_tile(color="black") +
        geom_tile(data=filter(X, name_q==name_t), color="black", size=1) +
        scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
        # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
        xlab( "CMap Target" )
}

composite_plot <- function(X) {
    ## DGE plot
    gg1 <- fplot( filter(X, source == "dge") ) +
        scale_x_discrete(position = "top") +
      # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(hjust=0, vjust=0.5),
              plot.margin = margin(r = 0.25, l = 0.25, unit = "in")) +
        ylab( "3' DGE Query" )

    ## L1000 plot
    gg2 <- fplot( filter(X, source == "l1000") ) +
        ylab( "L1000 Query" ) +
        # scale_fill_gradientn( colors=pal, name="Tau", limits=c(-100,100) ) +
        # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, name = "Tau", limits = c(-100, 100)) +
        theme(legend.position = "bottom",
              plot.margin = margin(r = 0.25, l = 0.25, unit = "in"))
#
#     cell_plot <- X %>%
#         distinct(drugQ, cell_id_query) %>%
#         # unchop(cell_id_query) %>%
#         mutate(
#             cell_id_query = map_chr(
#                 cell_id_query,
#                 ~if (is.null(.x) || length(.x) == 1) .x %||% "" else "multiple"
#             )
#         ) %>%
#         drop_na() %>%
#         ggplot(aes(drugQ, fill = cell_id_query)) +
#             geom_bar() +
#             coord_flip() +
#             theme_minimal() + theme_bold() +
#             scale_fill_brewer(palette = "Set2", name = "Query cell line") +
#             theme(
#                 # axis.ticks.y = element_blank(), axis.text.y = element_blank(),
#                 # axis.title.y = element_blank(), axis.title.x = element_blank(),
#                 axis.title = element_blank(), axis.ticks = element_blank(),
#                 axis.text.x = element_blank(), axis.text.y = element_blank()
#             )

    ## Summary plot
    S <- X %>% filter(name_q == name_t) %>%
        mutate_at("source", factor, levels=c("l1000","dge")) %>%
        mutate_at("source", fct_recode, `Self (3' DGE)`="dge", `Self (L1000)`="l1000")
    ggs <- ggplot( S, aes(x=name_t, y=source, fill=tau) ) +
        theme_minimal() + theme_bold() +
        geom_tile(color="black") + ylab("") +
        scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(-100,100) ) +
        # exec(scale_fill_gradientn, !!!cmap_nonlinear_colormap, guide = FALSE, limits = c(-100, 100)) +
        theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    ## Create the composite plot
    egg::ggarrange( gg1, ggs, gg2, heights=c(7.5,0.5,7.5), widths = c(7), draw=FALSE )
}

write_csv(
  R4, file.path(wd, "2b_data.csv")
)




ggcomp <- R4_completed %>%
    filter( source == "l1000" ) %>%
    group_nest( cutoff ) %>%
    mutate(
        data = map(
            data,
            bind_rows,
            filter( R4_completed, source == "dge")
        ) %>%
            map(composite_plot)
    )

ggsave(
  "fig2_summary_target_mean_query.pdf",
  ggcomp$data[[1]], width = 10, height = 18
)


 pwalk(
    ggcomp,
    function(cutoff, query_type, data, ...) {
        walk(
          file.path(
            wd,
            paste0("fig2_mcf7_target_", cutoff, c(".png", ".pdf"))
          ),
          ggsave, data, width=10, height=18 )
    }
)

# Test if self-similarity is significantly higher than other similarities

wilcox_self_similarities <- R3 %>%
    group_by(source, cutoff, lspci_id_q, name_q) %>%
    summarize(
        wilcox = wilcox.test(
            tau[lspci_id_q != lspci_id_t],
            mu = tau[lspci_id_q == lspci_id_t],
            alternative = "less"
        ) %>%
            list(),
        .groups = "drop"
    ) %>%
    mutate(
        p.value = map_dbl(wilcox, "p.value")
    )

wilcox_self_similarities %>%
  count(cutoff, source, significant = p.value < 0.05) %T>%
  write_csv(file.path(wd, "wilcox_table.csv"))

# 24 out of 32 DGE profiles show significant self-similarity

library(ggbeeswarm)

## Define transform that min-max normalizes data to between 0 and 1
## and then applies a logistic transformation
## Useful for plotting symmetrical data emphasizing points at the extreme ends
## of the distribution
logit_p_trans <- function(mi, ma) {
  scales::trans_new("logit_p",
                    transform = function (x) qlogis((x - mi) / (ma - mi)),
                    inverse   = function (x) ((ma - mi) * plogis(x)) + mi)
}

self_similarity_beeswarm_data <- R4 %>%
    filter(is.na(cutoff) | cutoff == 0.7) %>%
    mutate(
        self_similarity = if_else(
            lspci_id_q == lspci_id_t,
            "self_similarity",
            "cross_similarity"
        ) %>%
            fct_relevel("cross_similarity")
    ) %>%
    arrange(self_similarity)

self_similarity_beeswarm_plot <- self_similarity_beeswarm_data %>%
        ggplot(
            aes(tau, name_q, color = self_similarity)
        ) +
        geom_quasirandom(
            data = ~filter(.x, self_similarity == "cross_similarity"),
            groupOnX = FALSE
        ) +
        geom_point(
            data = ~filter(.x, self_similarity == "self_similarity")
        ) +
        facet_wrap(~source, nrow = 1) +
        scale_color_manual(
            values = c(
                self_similarity = "#FF0000",
                cross_similarity = "#00000088"
            ),
            guide = FALSE
        ) +
        scale_x_continuous(trans = logit_p_trans(-101, 101))

ggsave(
  file.path(wd, "similarity_beeswarm_separate.pdf"),
  self_similarity_beeswarm_plot, width = 6, height = 9
)


self_similarity_beeswarm_agg_plot <- self_similarity_beeswarm_data %>%
    ggplot(
        aes(source, tau, fill = self_similarity)
    ) +
    geom_quasirandom(
        data = ~filter(.x, self_similarity == "cross_similarity"),
        shape = 21,
        size = 1,
        color = "NA",
        method = "quasirandom",
        bandwidth = 0.2,
        width = 0.5
    ) +
    # geom_beeswarm(
    #     data = ~filter(.x, self_similarity == "cross_similarity"),
    #     shape = 21,
    #     color = "NA",
    #     priority = "random"
    # ) +
    scale_fill_manual(
        values = c(
            self_similarity = "#FF0000",
            cross_similarity = "#00000088"
        ),
        guide = FALSE
    ) +
    facet_wrap(~source, scales = "free", ncol = 1) +
    theme_light() +
    theme_bold() +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()
    ) +
    scale_y_continuous(position = "right", trans = logit_p_trans(-101, 101)) +
    labs(x = NULL, y = "Tau")


ggsave(
    file.path(wd, "fig2b_self_similarity_beeswarm_agg_logit.pdf"),
    self_similarity_beeswarm_agg_plot,
    width = 1.5, height = 7.5
)

ggsave(
  file.path(wd, "fig2b_self_similarity_beeswarm_agg_linear.pdf"),
  self_similarity_beeswarm_agg_plot +
    scale_y_continuous(position = "right"),
  width = 1.5, height = 7.5
)

self_similarity_stats <- R4 %>%
    group_by(source, lspci_id_q, name_q, cutoff) %>%
    arrange(desc(tau)) %>%
    summarize(
        rank = seq_len(n())[lspci_id_q == lspci_id_t],
        rank_percentile = rank / n(),
        .groups = "drop"
    )

wilcox_self_similarities_col_plot <- wilcox_self_similarities %>%
    mutate(log_p = -log10(p.value)) %>%
    ggplot(aes(x = log_p, y = name_q, fill = source)) +
    geom_col() +
    facet_wrap(~source, nrow = 1) +
    guides(fill = FALSE)

ggsave(
    file.path("self_similarity", "wilcox_cols.pdf"),
    wilcox_self_similarities_col_plot,
    width = 4, height = 8
)

wilcox_self_similarities_scatterplot <- wilcox_self_similarities %>%
    mutate(log_p = -log10(p.value)) %>%
    select(idQ, source, p.value) %>%
    spread(source, p.value) %>%
    ggplot(aes(x = DGE, y = L1000)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    coord_equal()



# Plot TAS vector similarity

# Use old version of TAS vectors that matches used lspci_ids
tas <- syn("syn20830942.4") %>%
    read_csv()

library(data.table)

lspci_id_name_map <- R2 %>%
    distinct(lspci_id = idQ, name = drugQ) %>%
    mutate(across(lspci_id, as.double))

tas_used <- tas %>%
    filter(fp_name == "morgan_normal", lspci_id %in% R2$idQ) %>%
    distinct(lspci_id, gene_id = entrez_gene_id, tas) %>%
    setDT()

tas_weighted_jaccard <- function(data_tas, query_id, min_n = 6) {
    query_tas <- data_tas[lspci_id == query_id, .(gene_id, tas = 11 - tas)]
    data_tas[
        ,
        .(lspci_id, gene_id, tas = 11 - tas)
    ][
        query_tas,
        on = "gene_id",
        nomatch = NULL
    ][
        ,
        mask := tas > 1 | i.tas > 1
    ][
        ,
        if (sum(mask) >= min_n) .(
            "tas_similarity" = sum(pmin(tas[mask], i.tas[mask])) / sum(pmax(tas[mask], i.tas[mask])),
            "n" = sum(mask),
            "n_prior" = .N
        ) else .(
            "tas_similarity" = numeric(),
            "n" = integer(),
            "n_prior" = integer()
        ),
        by = "lspci_id"
    ]
}

all_similarity <- tibble(lspci_id_1 = unique(tas_used$lspci_id)) %>%
    mutate(
        data = map(
            lspci_id_1,
            ~tas_weighted_jaccard(tas_used, .x) %>%
                rename(lspci_id_2 = lspci_id)
        )
    ) %>%
    unnest(data) %>%
    left_join(
        lspci_id_name_map %>%
            rename(name_1 = name, lspci_id_1 = lspci_id)
    ) %>%
    left_join(
        lspci_id_name_map %>%
            rename(name_2 = name, lspci_id_2 = lspci_id)
    ) %>%
    mutate(across(starts_with("name"), factor, levels = lvl)) %>%
    mutate(across(name_2, fct_rev))

tas_similarity_plot <- ggplot( all_similarity, aes(x=name_1, y=name_2, fill=tas_similarity) ) +
    theme_minimal() + theme_bold() +
    geom_tile(color="black") +
    geom_tile(data=filter(all_similarity, lspci_id_1==lspci_id_2), color="black", size=1) +
    # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(0, 1) ) +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
    # xlab( "CMap Target" ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(hjust=0, vjust=0.5),
          axis.title.y = element_blank())
    # ylab( "3' DGE Query")

dir.create("tas_similarity")
ggsave(
    file.path("tas_similarity", "tas_similarity_heatmap.pdf"),
    width = 8, height = 6
)

# Tas similarity Fig 6

fig6_drugs <- c(4763, 27399, 86140, 99986, 84200, 102661, 117785,
  101253, 36292, 57736, 100531, 14772, 97896, 76898,
  95012, 12104, 103943, 66433, 101674, 83449, 84593,
  96877, 52270, 72549, 96251, 82024, 99378, 66419,
  90309, 87501, 91047, 82566, 45745, 90255, 86536,
  99422, 96405, 75291, 92053, 97426, 20087, 66998,
  91759, 79027, 52760)

tas_used <- tas %>%
    filter(fp_name == "morgan_normal", lspci_id %in% fig6_drugs) %>%
    distinct(lspci_id, gene_id = entrez_gene_id, tas) %>%
    setDT()

all_similarity <- tibble(lspci_id_1 = unique(tas_used$lspci_id)) %>%
    mutate(
        data = map(
            lspci_id_1,
            ~tas_weighted_jaccard(tas_used, .x) %>%
                rename(lspci_id_2 = lspci_id)
        )
    ) %>%
    unnest(data) %>%
    mutate(across(starts_with("lspci_id"), factor, levels = fig6_drugs)) %>%
    mutate(across(lspci_id_2, fct_rev))

tas_similarity_plot <- ggplot( all_similarity, aes(x=lspci_id_1, y=lspci_id_2, fill=tas_similarity) ) +
    theme_minimal() + theme_bold() +
    geom_tile(color="black") +
    # geom_tile(data=filter(all_similarity, lspci_id_1==lspci_id_2), color="black", size=1) +
    # scale_fill_gradientn( colors=pal, guide=FALSE, limits=c(0, 1) ) +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90") +
    # xlab( "CMap Target" ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
# ylab( "3' DGE Query")

dir.create("tas_similarity")
ggsave(
    file.path("tas_similarity", "tas_similarity_fig6_heatmap.pdf"),
    tas_similarity_plot,
    width = 8, height = 6
)


## Query by cell-line


## Isolate the appropriate slice of data
## Aggregate across multiple entries to compute master similarity score
R2 <- gene_set_meta %>%
  filter(
    # cell_aggregate_method == "cells_aggregated" | is.na(cell_aggregate_method),
    # replicate_method == "replicates_aggregated" | is.na(replicate_method),
    # concentration_method == "concentration_aggregated" | is.na(concentration_method),
    is.na(stim)
  ) %>%
  rename(
    lspci_id_q = lspci_id, cells_q = cells
  ) %>%
  distinct() %>%
  inner_join(
    R %>%
      rename(cells_t = cell_id),
    by = c("gene_set_id" = "gene_set")
  ) %>%
  left_join(
    cmap_meta %>%
      distinct(pert_id, lspci_id_t = lspci_id),
    by = "pert_id"
  ) %>% {
    bind_rows(
      filter(., source == "l1000"),
      # Actually didn't perform cell aggregation for dge data... need to do here
      filter(., source == "dge") %>%
        group_by(source, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_t) %>%
        summarize(
          tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
            {.[order(abs(.))[2]]},
          .groups = "drop"
        ) %>%
        mutate(cell_aggregate_method = "cells_aggregated")
    )
  } %>%
  # For some reason, CMap sometimes returns multiple independent connectivities
  # of the same query and target compound. Probably replicate signatures
  # on their side?
  # Aggregating by taking the 33- or 66-percentile, whichever has
  # higher absolute value. Approach used by CMap to aggregate cell lines
  group_by(source, cell_aggregate_method, replicate_method, concentration_method, lspci_id_q, lspci_id_t, cells_q, cells_t, cutoff) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    .groups = "drop"
  ) %>%
  left_join(
    compound_names %>%
      rename(name_q = name),
    by = c("lspci_id_q" = "lspci_id")
  ) %>%
  left_join(
    compound_names %>%
      rename(name_t = name),
    by = c("lspci_id_t" = "lspci_id")
  )


R3 <- R %>% filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
    select( idQ, idT, tau, source, z_score_cutoff, query_type, cell_id_query ) %>%
    group_by( idQ, idT, source, z_score_cutoff, query_type, cell_id_query ) %>%
    summarize_at( "tau", ~.x[ which.max(abs(.x)) ] ) %>% ungroup() %>%
    mutate_at( c("idQ", "idT"), as.character ) %>%
    mutate_at( "source", toupper ) %>%
    mutate( drugT = dmap[idT], drugQ = dmap[idQ] )

R_torin <- R3 %>%
    filter(query_type == "per_cell_line") %>%
    mutate(drugQ = cell_id_query) %>%
    bind_rows(
        R2 %>%
            filter(query_type == "aggregated", source == "L1000", drugQ == "torin-1") %>%
            mutate(drugQ = "aggregated")
    ) %>%
    bind_rows(
        crossing(
            drugQ = "",
            drugT = unique(.[["drugT"]]),
            z_score_cutoff = unique(.[["z_score_cutoff"]])
        ) %>%
            mutate(tau = 0)
    )

plot_torin <- R_torin %>%
    group_nest(z_score_cutoff) %>%
    mutate(
        data = map(
            data,
            ~.x %>%
                mutate(
                    drugQ = fct_relevel(drugQ, "aggregated", ""),
                    border = case_when(
                        drugQ == "" ~ 0,
                        idQ == idT ~ 1,
                        TRUE ~ 0.4
                    ),
                    color = if_else(drugQ == "", NA_character_, "black")
                ) %>%
                ggplot(aes(x=drugT, y=drugQ, fill=tau, size=border) ) +
                scale_size_identity() +
                theme_minimal() + theme_bold() +
                geom_tile(aes(color = color)) +
                scale_color_identity() +
                # geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
                scale_fill_gradientn( colors=pal, limits=c(-100,100) ) +
                xlab( "Clue Target" ) +
                ylab("Query cell line")
        )
    )

pwalk(
    plot_torin,
    function(z_score_cutoff, data, ...) {
        walk(
            paste0("figSx_torin-1_per_cell_line_", z_score_cutoff, c(".pdf", ".png")),
            ggsave,
            data,
            width = 7, height = 6
        )
    }
)

## Target by cell-line

R_all <- file.path(pathData, "clue_results_combined.rds") %>% read_rds()

R4 <- R_all %>%
    filter( result_type == "pert", score_level == "cell" ) %>%
    pluck( "data", 1 ) %>%
    rename(idQ = lspci_id_query, idT = lspci_id_target, drugT = pert_iname)%>%
    filter(idT %in% names(dmap), idQ %in% names(dmap)) %>%
    select( idQ, idT, tau, source, z_score_cutoff, query_type, cell_id_query, cell_id_target ) %>%
    group_by( idQ, idT, source, z_score_cutoff, query_type, cell_id_target ) %>%
    summarize_at( "tau", ~.x[ which.max(abs(.x)) ] ) %>% ungroup() %>%
    mutate_at( c("idQ", "idT"), as.character ) %>%
    mutate_at( "source", toupper ) %>%
    mutate( drugT = dmap[idT], drugQ = dmap[idQ] )

R4_torin <- R4 %>%
    filter(query_type == "aggregated", idQ == 101253) %>%
    mutate(drugQ = str_to_upper(cell_id_target)) %>%
    bind_rows(
        R2 %>%
            filter(query_type == "aggregated", drugQ == "torin-1") %>%
            mutate(drugQ = "aggregated")
    ) %>%
    bind_rows(
        crossing(
            drugQ = "",
            drugT = unique(.[["drugT"]]),
            z_score_cutoff = unique(.[["z_score_cutoff"]]),
            source = unique(.[["source"]])
        ) %>%
            mutate(tau = 0)
    )

plot_torin2 <- R4_torin %>%
    group_nest(source, z_score_cutoff) %>%
    mutate(
        data = map(
            data,
            ~.x %>%
                mutate(
                    drugQ = fct_relevel(drugQ, "aggregated", ""),
                    border = case_when(
                        drugQ == "" ~ 0,
                        idQ == idT ~ 1,
                        TRUE ~ 0.4
                    ),
                    color = if_else(drugQ == "", NA_character_, "black")
                ) %>%
                ggplot(aes(x=drugT, y=drugQ, fill=tau, size=border) ) +
                scale_size_identity() +
                theme_minimal() + theme_bold() +
                geom_tile(aes(color = color)) +
                scale_color_identity() +
                # geom_tile(data=filter(X, idQ==idT), color="black", size=1) +
                scale_fill_gradientn( colors=pal, limits=c(-100,100) ) +
                xlab( "Clue Target" ) +
                ylab("Target cell line")
        )
    )

pwalk(
    plot_torin2,
    function(z_score_cutoff, data, source, ...) {
        walk(
            paste0("figSx_torin-1_target_cell_line_", source, "_", z_score_cutoff, c(".pdf", ".png")),
            ggsave,
            data,
            width = 7, height = 3.5
        )
    }
)


## Plot two matrix halfs against each other

halfs_data <- R4 %>%
  filter(name_q != name_t) %>%
  mutate(
    direction = if_else(lspci_id_q > lspci_id_t, "forward", "reverse"),
    across(c(name_q, name_t), as.character),
    names = map2(name_q, name_t, c) %>%
      map(sort),
    name_1 = map_chr(names, 1),
    name_2 = map_chr(names, 2),
    comp = names %>%
      map_chr(paste, collapse = "_")
  ) %>%
  select(source, cutoff, direction, name_1, name_2, comp, cells_q, tau) %>%
  pivot_wider(names_from = direction, values_from = c(cells_q, tau)) %>%
  mutate(
    same_cells = if_else(cells_q_forward == cells_q_reverse, "same", "different")
  )

outliers <- halfs_data %>%
  filter(is.na(cutoff) | cutoff == 0.7) %>%
  mutate(source = fct_recode(source, `DGE query` = "dge", `L1000 query` = "l1000")) %>%
  filter(
    abs(tau_reverse - tau_forward) > 180
  )

halfs_data %>%
  filter(is.na(cutoff) | cutoff == 0.7) %>%
  mutate(source = fct_recode(source, `DGE query` = "dge", `L1000 query` = "l1000")) %>%
  mutate(
    outlier = abs(tau_reverse - tau_forward) > 180
  ) %>%
  group_by(
    outlier
  ) %>%
  summarize(prop_same_cells = sum(same_cells == "same") / n())

library(GGally)

halfs_data %>%
  filter(is.na(cutoff) | cutoff == 0.7) %>%
  mutate(source = fct_recode(source, `DGE query` = "dge", `L1000 query` = "l1000")) %>%
  mutate(
    outlier = if_else(abs(tau_reverse - tau_forward) > 180, "disconcordant", "concordant") %>%
      as.factor()
  ) %>%
  ggplot(
    aes(outlier, fill = same_cells, by = outlier)
  ) +
  geom_bar(position = "fill") +
  geom_text(stat = "prop", position = position_fill(0.5))+
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  labs(fill = "Query cell types")

ggsave(file.path(wd, "concordant_forward_reverse_same_cell_bar.pdf"), width = 4, height = 4)

# # A tibble: 2 × 2
# outlier prop_same_cells
# <lgl>             <dbl>
#   1 FALSE             0.710
# 2 TRUE              0.303

library(ggrepel)
p <- ggplot(
  halfs_data %>%
    filter(is.na(cutoff) | cutoff == 0.7) %>%
    mutate(source = fct_recode(source, `DGE query` = "dge", `L1000 query` = "l1000")) %>%
    arrange(same_cells),
  aes(tau_forward, tau_reverse, color = same_cells)
) +
  geom_point(alpha = 0.8, shape = 16) +
  # facet_grid(vars(source), vars(same_cells)) +
  facet_wrap(~source) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  # geom_text_repel(
  #   aes(label = comp),
  #   data = ~.x %>%
  #     mutate(
  #       comp = if_else(abs(reverse - forward) > 180, comp, "")
  #     ),
  #     # filter(),
  #   max.overlaps = 1000
  # ) +
  labs(
    x = "Connectivity Query X Target Y", y = "Connectivity Query Y Target X",
    color = "Query cell types"
  )
  # lims(x = c(-130, 130), y = c(-130, 130))

ggsave(
  file.path(wd, "query_reciprocal_scatter.pdf"), p,
  width = 8, height = 4
)

# outliers <- halfs_data %>%
#   filter(
#     abs(reverse - forward) > 180
#   )
