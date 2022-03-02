library( tidyverse )
library( ggrepel )
library( ggbeeswarm )

pathData <- "~/data/DGEsig"


synapser::synLogin()
syn <- synExtra::synDownloader(pathData, ifcollision="overwrite.local")

compound_name_map <- syn("syn22035396.3") %>%
    read_rds() %>%
    filter(fp_name == "morgan_normal") %>%
    chuck("data", 1) %>%
    group_by(lspci_id) %>%
    slice(1L) %>%
    ungroup() %>%
    mutate(across(name, str_to_lower))

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
cndict <- c(a375 = "A375", a549 = "A549", ha1e = "HA1E", hcc515 = "HCC515",
            hepg2 = "HepG2", ht29 = "HT29", mcf7 = "MCF-7", pc3 = "PC3", vcap = "VCaP")

## Isolate the relevant data chunk
R0 <- syn("syn21907139.6") %>% read_rds()

R <- R0 %>%
    filter(result_type == "pert", score_level == "cell") %>%
    chuck("data", 1L) %>%
    filter( pert_type == "trt_cp", cell_id_query == "mcf7" ) %>%
    select(cellT = cell_id_target, drugT = pert_iname, tau, gene_set,
           idT = lspci_id_target, idQ = lspci_id_query) %>%
    filter( !is.na(idQ), !is.na(idT) ) %>%
    mutate_at( "cellT", recode, !!!cndict ) %>%
    left_join(
        compound_name_map %>%
            distinct(lspci_id, drugQ = name),
        by = c("idQ" = "lspci_id")
    )

## MCF7 @ Alvocidib / Flavopiridol / Palbo
S1 <- R %>% filter(
        idT == idQ,
        drugQ %in% c("alvocidib", "palbociclib")
    ) %>%
    inner_join( CM, by = c("cellT"="Short") ) %>%
    # mutate_at( c("drugQ","Label"), factor ) %>%
    rename( Tau = tau )
    # mutate(
    #     drugQ = if_else(
    #         drugQ == "alvocidib",
    #         paste0("alvocidib R", as.integer(as.factor(gene_set))),
    #         drugQ
    #     )
    # )

## Split by tau interval
SS1 <- S1 %>% mutate( region = fct_rev(cut(Tau, breaks=c(-100, 95, 100))) ) %>%
    split(.$region)

## Plotting elements
pal <- set_names( ggthemes::few_pal()(7), unique(S1$Tissue) )
etxt <- function(s, ...) {element_text( size = s, face = "bold", ... )}
theme_bold <- function() {
    theme(axis.text = etxt(12), axis.title = etxt(14),
          legend.text = etxt(12), legend.title = etxt(14))
}

## Additional parameters for each plot
fplot <- function( .df, isTop ) {
    gg <- ggplot(.df, aes(x=drugQ, y=Tau, color=Tissue, group=drugQ)) +
        theme_bw() + theme_bold() +
        geom_beeswarm(cex=3, beeswarmArgs=list(side=-1)) +
        geom_vline(xintercept=c(1.7,2.7), color="lightgray") +
        scale_color_manual( values=pal, drop=FALSE ) +
        scale_x_discrete( drop=FALSE, expand=expansion(add=c(0.3,0.6)),
                         name="3' DGE Query Signature" ) +
        geom_text_repel(aes(label=cellT), show.legend=FALSE,
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

ggs <- map2( SS1, c(TRUE,FALSE), fplot )
ggcomp <- egg::ggarrange( plots=ggs, ncol=1, heights=c(0.5,0.5), draw=FALSE )

ggsave("fig3.pdf", ggcomp, width=6, height=9)
ggsave("fig3.png", ggcomp, width=6, height=9)
