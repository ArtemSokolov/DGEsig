library( tidyverse )
library( seriation )   # For optimal leaf reordering
library(synapser)
library(here)
library(qs)
library(ggbeeswarm)
library(ggrepel)

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
R <- syn("syn26468923") %>%
  qread()

R_cells <- R %>%
  filter(result_type == "pert", score_level == "cell") %>%
  chuck("data", 1) %>%
  left_join(
    cmap_meta %>%
      distinct(lspci_id, pert_id)
  )

condition_conc_vars <- c("cells", "drug_id", "lspci_id", "stim", "stim_conc", "time")

M <- syn("syn25292310") %>%
  # syn("syn22000707") %>%
  qread() %>%
  unnest(meta)


cmap_gene_sets <- syn("syn25314203.4") %>%
  qread()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

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
cmap_returned <- filter(cmap_meta, pert_id %in% R_cells$pert_id)$lspci_id %>% unique()
dge_queried <- filter(dge_gene_sets, gene_set_id %in% R_cells$gene_set)$lspci_id %>% unique()
qcom <- intersect(cmap_returned, dge_queried) %>%
  na.omit()

## Identify overlap in cell line and perturbation between DGE query gene sets
## and targets returned by CMap


# overlapping_gene_sets <- R_cells %>%
#   semi_join(
#     gene_set_meta,
#     by = c(
#       "cell_id" = "cells",
#       "gene_set",
#       ""
#     )
#   )
overlapping_gene_sets <- gene_set_meta %>%
  filter(
    replicate_method == "replicates_aggregated",
    concentration_method == "concentration_aggregated"
  ) %>%
  semi_join(
    R_cells,
    by = c(
      "cells" = "cell_id",
      "gene_set_id" = "gene_set",
      "lspci_id"
    )
  ) %>%
  drop_na(lspci_id)

R_overlap <- R_cells %>%
  rename(cells_target = cell_id) %>%
  inner_join(
    overlapping_gene_sets %>%
      rename(cells_query = cells),
    by = c(
      "gene_set" = "gene_set_id",
      "lspci_id"
    )
  ) %>%
  group_by(
    lspci_id,
    cells_query,
    cells_target
  ) %>%
  summarize(
    tau = quantile(tau, c(0.67, 0.33), names = FALSE, na.rm = TRUE) %>%
      {.[order(abs(.))[2]]},
    .groups = "drop"
  )



## Cell line metadata
CM <- tribble(
  ~Short,    ~Tissue,    ~Subtype,
  "A375",   "Skin",     "melanoma",
  "A549",   "Lung",     "adeno",
  "HA1E",   "Kindey",   "normal",
  "HCC515", "Lung",     "adeno",
  "HepG2",  "Liver",    "hcc",
  "HT29",   "Colon",    "adeno",
  "MCF-7",  "Breast",   "lumA",
  "PC3",    "Prostate", "small cell",
  "VCaP",   "Prostate", "AR-V7" ) %>%
  mutate( Label = glue::glue("{Short} ({Tissue}/{Subtype})") )

## Dictionary for cell name standardization
cndict <- c(A375 = "A375", A549 = "A549", HA1E = "HA1E", HCC515 = "HCC515",
            HEPG2 = "HepG2", HT29 = "HT29", MCF7 = "MCF-7", PC3 = "PC3", VCAP = "VCaP")

R_plotting <- R_overlap %>%
  mutate(across(c(cells_query, cells_target), recode, !!!cndict)) %>%
  arrange(tau) %>%
  group_by(
    lspci_id,
    cells_query
  ) %>%
  mutate(
    rank = factor(rev(seq_len(n())))
  ) %>%
  ungroup() %>%
  left_join(
    compound_names
  ) %>%
  mutate(
    region = fct_rev(
      cut(tau, breaks = c(-100, 95, 100))
    )
  ) %>%
  left_join(
    CM,
    by = c("cells_target" = "Short")
  )

## Plotting elements
pal <- set_names( ggthemes::few_pal()(7), unique(CM$Tissue) )
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
  theme(axis.text = etxt(12), axis.title = etxt(14),
        legend.text = etxt(12), legend.title = etxt(14))
}

## Additional parameters for each plot
fplot <- function( .df, isTop ) {
  gg <- ggplot(.df, aes(x=name, y=tau, color=Tissue, group=name)) +
    theme_bw() + theme_bold() +
    geom_beeswarm(cex=3, beeswarmArgs=list(side=-1)) +
    scale_color_manual( values=pal, drop=FALSE ) +
    scale_x_discrete( drop=FALSE, expand=expansion(add=c(0.3,0.6)),
                      name="3' DGE Query Signature" ) +
    geom_text_repel(aes(label=cells_target), show.legend=FALSE,
                    fontface="bold", nudge_x = 0.3,
                    direction="y", segment.color=NA ) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

  if( isTop == TRUE )
    gg <- gg + theme(axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.x = element_blank(),
                     legend.position="top")
  else
    gg <- gg + theme(axis.text.x = element_text(hjust=0.3)) +
      guides( color=FALSE )
  gg
}

ggs <- map2( split(R_plotting, R_plotting$region), c(TRUE,FALSE), fplot )
ggcomp <- egg::ggarrange( plots=ggs, ncol=1, heights=c(0.5,0.5), draw=FALSE )

ggsave("fig3.pdf", ggcomp, width=6, height=9)
ggsave("fig3.png", ggcomp, width=6, height=9)

p <- ggplot(R_plotting, aes(name, cells_target, fill = tau)) +
  geom_raster() +
  scale_fill_distiller(limits = c(-100, 100), palette = "RdBu", trans = logit_p_trans(-101, 101))
ggsave(
  "fig3c_hm.pdf", p, width = 15, height = 5
)

p <- ggplot(R_plotting, aes(name, rank, fill = Tissue)) +
  geom_raster() +
  geom_text(aes(label = cells_target)) +
  scale_fill_manual( values=pal, drop=FALSE )

p <- ggplot(R_plotting, aes(rank, name, fill = cells_target == "MCF-7")) +
  geom_tile(color = "black") +
  geom_text(aes(label = cells_target)) +
  scale_fill_manual( values=c(`TRUE` = "red", `FALSE` = "white") ) +
  guides(fill = "none") +
  theme_minimal() +
  theme_bold() +
  labs(x = "Tau rank", y = "Compound")
ggsave(
  "fig3c_rank.pdf", p, width = 8, height = 5
)


# library(scales)

# asinh_breaks <- function(r) {
#   lmin <- round(log10(abs(r[1]))) * sign(r[1])
#   lmax <- round(log10(r[2]))
#   lbreaks <- seq(lmin, lmax, by = 1)
#   breaks <- 10 ^ lbreaks
# }
#
# asinh_trans <- function() {
#   trans_new("asinh",
#             transform = asinh,
#             inverse   = sinh,
#             breaks = asinh_breaks)
# }
#
# sinh_trans <- function(div = 1) {
#   trans_new("asinh",
#             transform = function (x) sinh(x / 100),
#             inverse   = function (x) 100 * asinh(x))
#             # breaks = asinh_breaks)
# }

## Define transform that min-max normalizes data to between 0 and 1
## and then applies a logistic transformation
## Useful for plotting symmetrical data emphasizing points at the extreme ends
## of the distribution
logit_p_trans <- function(mi, ma) {
  scales::trans_new("logit_p",
            transform = function (x) qlogis((x - mi) / (ma - mi)),
            inverse   = function (x) ((ma - mi) * plogis(x)) + mi)
}

# MCF-7 vs all others t-test
r <- R_plotting %>%
  group_by(name) %>%
  summarize(
    res = t.test(
      tau[cells_target != "MCF-7"], mu = tau[cells_target == "MCF-7"],
      alternative = "less"
    ) %>%
      broom::tidy() %>%
      list()
  ) %>%
  unnest(res)

# ggplot(R_plotting, aes(x = tau/2, y = name, color = Tissue)) +
p <- ggplot(R_plotting, aes(x = tau, y = name, color = if_else(cells_target == "MCF-7", "MCF-7", "other cell lines"))) +
  geom_quasirandom(groupOnX = FALSE, width = 0.2, data = ~filter(.x, cells_target != "MCF-7")) +
  geom_point(data = ~filter(.x, cells_target == "MCF-7"), size = 2) +
  # geom_text(
  #   aes(x = -100, y = name, label = signif(p.value, 2), fontface = if_else(p.value < 0.05, "bold", "plain")),
  #   data = r,
  #   hjust = 0,
  #   inherit.aes = FALSE
  # ) +
  scale_color_manual(
    values = c(`other cell lines` = "black", `MCF-7` = "red")
  ) +
  # guides(color = "none") +
  scale_x_continuous(limits = c(-100, 100), trans = logit_p_trans(-101, 101)) +
  facet_grid(rows = vars(name), scales = "free") +
  theme_bw() +
  theme_bold() +
  theme(
    strip.text = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0, "pt"),
    plot.margin = unit(c(8, 0, 8, 8), "pt"),
    legend.position = "top"
  ) +
  labs(y = "3' DGE Query Signature", x = "Tau", color = NULL)
  # egg::ggarrange(
  #   geom_text(
  #     aes(x = -100, y = name, label = signif(p.value, 2), fontface = if_else(p.value < 0.05, "bold", "plain")),
  #     data = r,
  #     hjust = 0,
  #     inherit.aes = FALSE
  #   ) +
  # )
  # scale_x_continuous(trans = sinh_trans(10000))
  # geom_text_repel(
  #   aes(label = cells_target), data = ~filter(.x, cells_target == "MCF-7")
  # )

pval_plot <- ggplot(r) +
  geom_text(
    aes(x = 0, y = name, label = signif(p.value, 2), fontface = if_else(p.value < 0.05, "bold", "plain")),
    hjust = 0,
    inherit.aes = FALSE
  ) +
  facet_grid(rows = vars(name), scales = "free") +
  # theme_bw() +
  theme_bold() +
  theme_void() +
  theme(
    panel.spacing = unit(0, "pt"),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 8), "pt")
  ) +
  coord_cartesian(xlim=c(0, 0.01))


ggcomp <- egg::ggarrange( p, pval_plot, ncol=2, widths = c(8, 3.5), draw=FALSE )

1 / mean(1 / r$p.value)

ggsave("fig3c.pdf", ggcomp, width = 5, height = 6)
