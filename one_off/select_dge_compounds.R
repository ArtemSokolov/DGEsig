# Selecting compounds from the LINCS1 screening library that were assayed
# with L1000 and a bunch of other assays

library(tidyverse)
library(data.table)
library(fst)
library(here)
library(cmapR)
library(dtplyr)
library(qs)
library(cluequery)

synapser::synLogin()

syn <- synExtra::synDownloader(here("data"))

lincs1_raw <- syn("syn25582054") %>%
  read_csv()

perturbation_meta <- syn("syn21547097") %>%
  read_csv()

signature_meta <- syn("syn21547101") %>%
  fread()

cmap_gene_meta <- syn("syn21547102") %>%
  fread()

lspci_id_vendor_id_map <- syn("syn24183248") %>%
  fread()

dge_meta <- syn("syn25292310") %>%
  qread() %>%
  unnest(meta)

niepel_compound_meta <- syn("syn25586657") %>%
  read_tsv()

niepel_response_classes <- syn("syn25586607") %>%
  read_tsv() %>%
  left_join(
    niepel_compound_meta %>%
      distinct(`Drug Name` = DrugName, FacilityID = HMSLid),
    by = "Drug Name"
  )

# In order to figure out which compounds are succesfully returned as targets
# by Clue, checking results of a recent (03/05/21) query of Clue

clue_res <- clue_query_download("604237ff198678001150fe96") %>%
  clue_parse_result(score_level = "cell", score_type = "tau", result_type = "pert")

clue_res_lspci_id <- clue_res %>%
  inner_join(
    distinct(perturbation_meta, lspci_id, pert_id),
    by = "pert_id"
  )

lincs1_compounds <- lincs1_raw %>%
  left_join(
    lspci_id_vendor_id_map %>%
      select(lspci_id, vendor_id) %>%
      as_tibble(),
    by = c("FacilityID" = "vendor_id")
  ) %>%
  mutate(
    cmap_mcf7 = lspci_id %in% {
      signature_meta %>%
        filter(cell_id == "MCF7") %>%
        pull(lspci_id)
    },
    clue_mcf7_return = lspci_id %in% {
      clue_res_lspci_id %>%
        filter(cell_id == "MCF7") %>%
        pull(lspci_id)
    }
  ) %>%
  left_join(
    niepel_response_classes %>%
      filter(`Cell line` == "MCF7") %>%
      distinct(FacilityID, Concentration, `Response class`, `GR value (72h)`) %>%
      arrange(FacilityID, Concentration) %>%
      group_by(FacilityID) %>%
      summarize(
        response_class = paste(
          paste(Concentration, `Response class`, sep = ":"),
          collapse = " "
        ),
        gr_value_72h = paste(
          paste(Concentration, signif(`GR value (72h)`, digits = 2), sep = ":"),
          collapse = " "
        ),
        .groups = "drop"
      ),
    by = "FacilityID"
  )

lincs1_compounds_cmap <- lincs1_compounds %>%
  inner_join(
    signature_meta %>%
      distinct(sig_id, lspci_id) %>%
      as_tibble(),
    by = "lspci_id"
  )

cmap_paths <- list(
  GSE92742 = "syn21551046",
  GSE70138 = "syn21551043",
  LINCS_2020 = "syn25050283"
) %>%
  enframe("dataset", "synid") %>%
  mutate(
    path = map_chr(synid, syn),
    rids = map(
      path,
      read_gctx_ids, dim = "row"
    ),
    cids = map(
      path,
      read_gctx_ids, dim = "col"
    )
  )


cmap_mats <- cmap_paths %>%
  mutate(
    mat = pmap(
      list(path, rids, cids),
      ~parse_gctx(..1, cid = which(..3 %in% lincs1_compounds_cmap[["sig_id"]]))
    )
  )

cmap_mat <- cmap_mats %>%
  pull(mat) %>%
  map(~.x@mat[order(rownames(.x@mat)), ]) %>%
  {do.call(cbind, .)} %>%
  # Some signatures duplicated, removing them
  {.[, !duplicated(colnames(.))]}

cmap_df <- cmap_mat %>%
  as_tibble(rownames = "pr_gene_id") %>%
  mutate(across(pr_gene_id, as.numeric)) %>%
  inner_join(
    cmap_gene_meta %>%
      distinct(pr_gene_id, entrez_id),
    by = "pr_gene_id"
  ) %>%
  drop_na(entrez_id)

dir.create(here("lincs1"))
write_fst(
  cmap_df,
  here("lincs1", "cmap_signatures_lincs1.fst")
)
# cmap_df <- read_fst(here("lincs1", "cmap_signatures_lincs1.fst"), as.data.table = TRUE)

# Selecting signatures closest to 24h

lincs1_chosen_signature_meta <- signature_meta %>%
  as_tibble() %>%
  filter(sig_id %in% colnames(cmap_df), cell_id == "MCF7") %>%
  mutate(time_dist = abs(24 - pert_time), pert_dose = signif(pert_dose, digits = 2)) %>%
  arrange(lspci_id, pert_dose) %>%
  group_by(lspci_id) %>%
  filter(time_dist == min(time_dist)) %>%
  mutate(
    cmap_doses = paste(unique(pert_dose), collapse = " ")
  ) %>%
  filter(pert_dose == max(pert_dose)) %>%
  ungroup()

lincs1_chosen_signatures <- lazy_dt(cmap_df) %>%
  select(pr_gene_id, any_of(lincs1_chosen_signature_meta$sig_id)) %>%
  inner_join(
    distinct(cmap_gene_meta, pr_gene_id, entrez_id),
    by = "pr_gene_id"
  ) %>%
  select(-pr_gene_id) %>%
  as.data.table() %>%
  melt(id.vars = "entrez_id", variable.name = "sig_id", value.name = "zscore") %>%
  merge(
    distinct(lincs1_chosen_signature_meta, sig_id, lspci_id, cell_id),
    by = "sig_id", all = FALSE, allow.cartesian = TRUE
  ) %>% {
    .[
      # Aggregate across replicates
      ,
      .(zscore = mean(zscore)),
      keyby = c("entrez_id", "lspci_id", "cell_id")
    ]
    # ][
    #   # Then across cell lines
    #   ,
    #   .(
    #     zscore = quantile(zscore, c(0.67, 0.33), names = FALSE) %>%
    #       {.[order(abs(.))[2]]}
    #   ),
    #   keyby = c("entrez_id", "lspci_id")
    # ]
  }

lincs1_chosen_signatures_diff <- lazy_dt(lincs1_chosen_signatures) %>%
  filter(abs(zscore) > qnorm(0.95)) %>%
  mutate(
    de = fifelse(zscore > 0, "up", "down")
  ) %>%
  as.data.table()

lincs1_compound_ranking <- lincs1_compounds %>%
  left_join(
    lincs1_chosen_signatures_diff %>%
      as_tibble() %>%
      dplyr::count(lspci_id, name = "n_deg_cmap"),
    by = "lspci_id"
  ) %>%
  left_join(
    lincs1_chosen_signature_meta %>%
      as_tibble() %>%
      distinct(lspci_id, cmap_doses),
    by = "lspci_id"
  ) %>%
  arrange(desc(n_deg_cmap)) %>%
  mutate(
    already_dge_profiled_mcf7 = lspci_id %in% {
      dge_meta %>%
        filter(cells == "MCF7") %>%
        pull(lspci_id)
    }
  ) %>%
  relocate(n_deg_cmap, already_dge_profiled_mcf7, cmap_doses, .after = clue_mcf7_return)

fwrite(
  lincs1_compound_ranking, here("lincs1_priorities.csv")
)
