library( tidyverse )
library( seriation )   # For optimal leaf reordering
library( cowplot )
library(synExtra)
library(data.table)
library(qs)

synapser::synLogin()

syn <- synDownloader("~/data", .cache = TRUE)

compound_names <- syn("syn26260344") %>%
  read_csv() %>%
  select(lspci_id, name) %>%
  drop_na() %>%
  group_by(lspci_id) %>%
  slice(1) %>%
  ungroup()

dge_gene_sets <- syn("syn25303778") %>%
  qread()

GS0 <- syn("syn22105667") %>%
  read_csv()

GS <- GS0 %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    is.na(stim),
    str_detect(coll(cell_id, ignore_case = TRUE), "rencell")
  ) %>%
    group_by( DrugID=lspci_id ) %>% summarize( Set = list(entrezgene_id) )

R0 <- syn("syn26468923") %>%
  qread()

## Perform additional wrangling
relevant_gene_sets <- dge_gene_sets %>%
  filter(
    concentration_method == "concentration_aggregated",
    replicate_method == "replicates_aggregated",
    is.na(stim),
    str_detect(coll(cells, ignore_case = TRUE), "rencell")
    # pert_type == "trt_cp"
  )

R <- R0 %>%
  filter(result_type == "pert", score_level == "summary") %>%
  chuck("data", 1) %>%
  inner_join(
    relevant_gene_sets %>%
      select(idQ = lspci_id, gene_set = gene_set_id),
    by = "gene_set"
  ) %>%
  select( idQ, idT = pert_id, tau )

# For GSEA
expr <- syn("syn25303743") %>%
  fread()

expr_agg <- expr[
  replicate_method == "replicates_aggregated" &
    concentration_method == "concentration_aggregated" &
    stim == "" &
    str_detect(coll(cells, ignore_case = TRUE), "rencell")
]

library(fgsea)
library(msigdbr)

# all_msigdbr <- msigdbr(category = "H") %>%
#   bind_rows(msigdbr(category = "C5", subcategory = "BP"))
all_msigdbr <- msigdbr() %>%
  filter(
    str_starts(gs_subcat, fixed("CP")) |
      gs_cat == "H"
  )

hallmark_gene_sets <- all_msigdbr %>%
  group_nest(gs_name) %>% {
    set_names(
      map(.$data, "human_ensembl_gene"),
      .$gs_name
    )
  }

# expr_for_gsea <- expr_agg %>%
#   filter(is.finite(log2FoldChange)) %>%
#   group_nest(lspci_id)

# fgsea_res <- expr_for_gsea %>%
  # head(n = 1) %>%
  # unnest(data) %>%
fgsea_res <- expr_agg %>%
  filter(is.finite(log2FoldChange)) %>%
  group_by(lspci_id) %>%
  summarize(
    # gsea_res = map(
    #   gene_set_table,
    #   ~fgseaMultilevel(
    #     hallmark_gene_sets, set_names(.x$log2FoldChange, .x$ensembl_gene_id)
    #   )
    # )
    res = {
      message(lspci_id[1])
      fgseaMultilevel(
        hallmark_gene_sets, set_names(log2FoldChange, ensembl_gene_id),
        nproc = 2
      ) %>%
        list()
    },
    .groups = "drop"
  )

library(qs)
qsave(
  fgsea_res,
  "fgsea_res2.qs"
)
fgsea_res <- qread("fgsea_res2.qs")

library(VIM)

fgsea_res_mat <- reduce2(
  fgsea_res$res, fgsea_res$lspci_id,
  function(agg, res, lspci_id) {
    message(lspci_id)
    lspci_id_s <- sym(as.character(lspci_id))
    full_join(
      agg,
      # select(res, pathway, {{lspci_id_s}} := padj),
      select(res, pathway, {{lspci_id_s}} := NES),
      by = "pathway"
    )
  },
  .init = tibble(pathway = character())
) %>%
  column_to_rownames("pathway") %>%
  as.matrix() %>%
  kNN( imp_var = FALSE, trace = TRUE)
  # kNN(numFun = laeken::weightedMean, weightDist = TRUE, imp_var = FALSE, trace = TRUE)

library(pheatmap)
withr::with_pdf(
  "fgsea_res_padj.pdf",
  pheatmap(
    fgsea_res_mat %>%
      as.matrix() %>% {
        -log10(.)
      }
  )
)
withr::with_pdf(
  "fgsea_res_nes.pdf",
  pheatmap(
    fgsea_res_mat %>%
      as.matrix()
  )
)

fgsea_cor <- cor(fgsea_res_mat, method = "pearson", use = "all.obs") %>%
  {
    na_mat <- is.na(.)
    .[-which(rowSums(na_mat) > 1), -which(colSums(na_mat) > 1)]
  }

fgsea_cor_df <- fgsea_cor %>%
  as_tibble(rownames = "DrugID1") %>%
  pivot_longer(-DrugID1, names_to = "DrugID2", values_to = "Correlation") %>%
  mutate(across(starts_with("DrugID"), as.numeric))

## Compose tau profiles for each drug
## Join against gene sets
V <- R %>% group_by( DrugID=idQ ) %>% arrange(idT) %>%
    summarize( Vals = list(set_names(tau, idT)) ) %>%
    inner_join(GS, by="DrugID")

## Ensure that all tau values are in the same order
stopifnot( map(V$Vals, names) %>% map_lgl(identical, .[[1]]) %>% all )

## Compute pair-wise similarity for all query drugs
tausim <- function( v1, v2 ) cor(v1,v2,method="pearson", use="complete.obs")
jcrdsim <- function(gs1, gs2)
    length(intersect(gs1, gs2)) / length(union(gs1, gs2))
SM <- crossing(rename_all(V, str_c, "1"),
               rename_all(V, str_c, "2")) %>%
    mutate(TauSim  = map2_dbl(Vals1, Vals2, tausim),
           JcrdSim = map2_dbl(Set1, Set2, jcrdsim))

## Perform hierarchical clustering with optimal leaf reordering
## Fixes row order based on the clustering results
## .df - data frame, .sim - similarity column
simorder <- function( .df, .sim )
{
    ## Compose the distance matrix
    DM <- .df %>% select( DrugID1, DrugID2, {{.sim}} ) %>%
        spread( DrugID2, {{.sim}} ) %>% as.data.frame() %>%
        column_to_rownames("DrugID1") %>% dist()

    ## Perform hierarchical clustering with optimal leaf reordering
    lvl <- hclust(DM) %>% reorder(DM) %>%
        dendextend::order.hclust() %>% labels(DM)[.]

    ## Fix order of rows and column based on clustering results
    .df  %>%
      left_join(
        compound_names %>%
          select(DrugID1 = lspci_id, DrugName1 = name),
        by = "DrugID1"
      ) %>%
      left_join(
        compound_names %>%
          select(DrugID2 = lspci_id, DrugName2 = name),
        by = "DrugID2"
      ) %>%
      transmute(DrugID1 = factor(DrugID1, lvl),
                      DrugID2 = factor(DrugID2, rev(lvl)),
                      Similarity = {{.sim}}, DrugName1, DrugName2)
}

## Plotting elements
pal  <- rev(RColorBrewer::brewer.pal(n=7, name="RdBu"))
##palj <- gray.colors(7, start=1, end=0)
palj <- pal[4:7]
ebl <- element_blank
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}

## Plots a similarity matrix
simplot <- function( .df )
{
    ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
        theme_minimal() +
        geom_tile(color="gray") +
##        geom_rect(aes(xmin=8.5, xmax=13.5, ymin=ndg-8.5+1, ymax=ndg-13.5+1),
##                  fill=NA, color="black", size=1) +
        theme(
          axis.text = ebl(),
              axis.title=ebl(),
              legend.text = etxt(12), legend.title=etxt(14),
              plot.margin = unit(c(0.5,0.5,0.5,2), "cm"))
}

## Plots a zoom panel
zoomplot <- function( .df )
{
    ggplot( .df, aes(DrugID1, DrugID2, fill=Similarity) ) +
        theme_minimal() + geom_tile(color="gray") +
        scale_y_discrete( labels = function(x) dnm[x] ) +
        theme(axis.title=ebl(), axis.text.x=ebl(), axis.ticks.x=ebl(),
              axis.text.y=etxt(12), plot.background=element_rect(color="black", size=2),
              panel.grid.minor=ebl(), panel.grid.major=ebl())
}

## Compose similarity matrices
XTau <- simorder( SM, TauSim )
XJcrd <- simorder( SM, JcrdSim ) %>%
    mutate(across(DrugID1, fct_relevel, levels(XTau$DrugID1)),
           across(DrugID2, fct_relevel, levels(XTau$DrugID2)))
ndg <- length(unique(SM$DrugID2))

XCor <- simorder(fgsea_cor_df, Correlation)
XCor_ordered <- XCor %>%
  mutate(across(DrugID1, fct_relevel, levels(XJcrd$DrugID1)),
         across(DrugID2, fct_relevel, levels(XJcrd$DrugID2)))

## Plot zoom facets
dnm <- set_names(c("bortezomib", "mg132", "fedratinib", "staurosp. agl.", "nilotinib"),
                 c("57736", "36292", "97896", "14772", "100531"))
ZTau <- XTau %>% filter( DrugID1 %in% names(dnm), DrugID2 %in% names(dnm) )
ZJcrd <- XJcrd %>% filter( DrugID1 %in% names(dnm), DrugID2 %in% names(dnm) )

## Plot similarity matrices
ggjcrd <- simplot(XJcrd) + scale_fill_gradientn( colors=palj, limits=c(0,1), name="Jaccard\nSimilarity" )
ggtau  <- simplot(XTau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Pearson\nCorrelation" )
ggcor  <- simplot(XCor)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Pearson\nCorrelation" )
ggcor_ordered  <- simplot(XCor_ordered)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), name="Pearson\nCorrelation" )

ggsave(
  "fig6_fgsea_cor_nes.pdf", ggcor, width = 8, height = 6
)
ggsave(
  "fig6_fgsea_cor_nes_ordered.pdf", ggcor_ordered, width = 8, height = 6
)

XCor_vs_XTau <- XCor %>%
  rename(Correlation = Similarity) %>%
  inner_join(
    XTau, by = c("DrugID1", "DrugID2", "DrugName1", "DrugName2")
  )

ggsave(
  "XCor_vs_XTau.pdf",
  ggplot(XCor_vs_XTau, aes(Similarity, Correlation)) +
    geom_point()
)

cor.test(
  XCor_vs_XTau$Correlation, XCor_vs_XTau$Similarity
)
# t = 82.646, df = 3247, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.8118793 0.8340570
# sample estimates:
#   cor
# 0.8232821

## Plot zoom facets
gzjcrd <- zoomplot(ZJcrd) + scale_fill_gradientn( colors=palj, limits=c(0,1), guide=FALSE )
gztau  <- zoomplot(ZTau)  + scale_fill_gradientn( colors=pal, limits=c(-1,1), guide=FALSE )

## Put everything together
gg <- withr::with_pdf("test.pdf", egg::ggarrange( plots=list(ggjcrd, ggtau), ncol=2,
                     labels=c(" A"," B"), padding=unit(2,"line"),
                     label.args = list(gp = grid::gpar(font = 4, cex = 4)) )) #%>%
##    ggdraw() %>%
##    + draw_plot( gztau, .68, .7, .17, .25 ) %>%
##    + draw_plot( gzjcrd, .18, .7, .17, .25 )
ggsave( "fig6.pdf", gg, width=17, height=7 )
ggsave( "fig6.png", gg, width=17, height=7 )

clust1 <- compound_names %>%
  filter(lspci_id %in% c(102509, 84708, 87368, 97675))

clusts <- tribble(
  ~"clust", ~"index",
  "1", 5:8,
  "2", 12:18
) %>%
  unchop(index) %>%
  mutate(
    lspci_id = levels(XTau$DrugID1)[index] %>%
      as.double()
  ) %>%
  left_join(
    compound_names
  )

library(data.table)
tas <- syn("syn26260405") %>%
  read_csv() %>%
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

tas_used <- tas %>%
  filter(
    lspci_id %in% as.double(levels(XTau$DrugID1))
  )

all_tas_similarity <- tibble(lspci_id_1 = as.double(levels(XTau$DrugID1))) %>%
  mutate(
    data = map(
      lspci_id_1,
      ~tas_weighted_jaccard(tas_used, .x) %>%
        rename(lspci_id_2 = lspci_id)
    )
  ) %>%
  unnest(data) %>%
  left_join(
    compound_names %>%
      rename(name_1 = name, lspci_id_1 = lspci_id)
  ) %>%
  left_join(
    compound_names %>%
      rename(name_2 = name, lspci_id_2 = lspci_id)
  ) %>%
  mutate(
    DrugID1 = factor(as.character(lspci_id_1), levels = levels(XTau$DrugID1)),
    DrugID2 = factor(as.character(lspci_id_2), levels = rev(levels(XTau$DrugID1))),
    Similarity = tas_similarity
  )

ggtas  <- simplot(all_tas_similarity)  +
  scale_fill_viridis_c( limits=c(0,1), name="Target\nJaccard\nSimilarity" ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE)

ggsave( "fig6C.png", ggtas, width=9, height=7 )

ggtau2 <- ggtau +
  geom_text(
    aes(label = text),
    data = all_tas_similarity %>%
      filter(tas_similarity > 0.5) %>%
      mutate(
        text = "*"
      )
  )

## Put everything together
gg <- egg::ggarrange( plots=list(ggjcrd, ggtau2, ggtas), ncol=2,
                      labels=c(" A"," B", " C"), padding=unit(2,"line"),
                      label.args = list(gp = grid::gpar(font = 4, cex = 4)) ) #%>%
##    ggdraw() %>%
##    + draw_plot( gztau, .68, .7, .17, .25 ) %>%
##    + draw_plot( gzjcrd, .18, .7, .17, .25 )
ggsave( "fig6.pdf", gg, width=17, height=15 )
ggsave( "fig6.png", gg, width=17, height=15 )

tau_vs_tas <- XTau %>%
  inner_join(
    all_tas_similarity %>%
      rename(TASSimilarity = Similarity),
    by = c("DrugID1", "DrugID2")
  )

p <- tau_vs_tas %>%
  ggplot(aes(Similarity, TASSimilarity)) +
  geom_point()

ggsave(
  "tau_vs_tas.pdf",
  p, width = 5, height = 5
)
