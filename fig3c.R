library(synExtra)
library(tidyverse)
library(here)
library(seriation)
library(pheatmap)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data/DGEsig")

wd <- here("fig3")
dir.create(wd, showWarnings = FALSE)

theme_set(theme_bw())

compound_name_map <- syn("syn22035396.3") %>%
  read_rds() %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1)

cmap_gene_meta <- syn("syn21547102") %>%
  read_csv()

clue_res_combined <- syn("syn21907166.4") %>%
  read_rds()
# clue_res_combined_2 <- clue_res_combined

cdk_lspci_ids <- compound_name_map %>%
  filter(name %in% c("Palbociclib", "Ribociclib", "Alvocidib", "Abemaciclib")) %>%
  pull(lspci_id) %>%
  unique()

cdk_results <- clue_res_combined %>%
  filter(result_type == "pert", score_level == "summary") %>%
  chuck("data", 1L) %>%
  filter(
    source == "dge", lspci_id_query %in% cdk_lspci_ids
    )

cdk_result_mat <- cdk_results %>%
  select(pert_id, name_query, tau) %>%
  group_by(name_query, pert_id) %>%
  summarize(
    tau = quantile(tau, c(0.33, 0.66), na.rm = TRUE, names = FALSE) %>% {
      .[which.max(abs(.))]
    },
    .groups = "drop"
  ) %>%
  pivot_wider(pert_id, names_from = name_query, values_from = tau) %>%
  drop_na() %>%
  column_to_rownames("pert_id") %>%
  as.matrix()

cdk_result_cor <- cdk_result_mat %>%
  cor()

cdk_result_clust <- cdk_result_cor %>%
  dist() %>%
  hclust() %>%
  reorder(dist(cdk_result_cor))

cdk_result_cor_heatmap <- pheatmap(
  cdk_result_cor,
  cluster_cols = cdk_result_clust,
  cluster_rows = cdk_result_clust,
  color = RColorBrewer::brewer.pal(5, "Reds") %>%
    colorRampPalette() %>% {.(100)}
)

cdk_result_cor_heatmap <- ggplot(
  cdk_result_cor %>%
    as_tibble(rownames = "query_x") %>%
    pivot_longer(-query_x, names_to = "query_y", values_to = "correlation") %>%
    mutate(
      across(
        starts_with("query"),
        factor,
        levels = cdk_result_clust$labels
      )
    ),
  aes(query_x, query_y, fill = correlation)
) +
  geom_tile() +
  scale_fill_distiller(palette = "Reds", direction = 1) +
  labs(x = NULL, y = NULL, fill = "CMap result\ncorrelation") +
  theme_minimal()

ggsave(
  file.path(wd, "fig3c_cdk_inhib_cmap_correlation.pdf"),
  cdk_result_cor_heatmap$gtable,
  width = 3.5, height = 3
)
