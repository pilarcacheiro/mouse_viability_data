### retrieve and curate single gene knockout viability annotations
### from the IMPC (https://www.mousephenotype.org/) and
### MGI (https://www.informatics.jax.org/) resources


# libraries and auxiliary scripts and files -------------------------------

if (!require("tidyverse")) install.packages("tidyverse")
library("tidyverse")

load("lethal_terms.RData")

### this is a vector containing a set of mammalian phenotypes ids
### extracted from Dickinson et al.(PMID: 27626380): embryonic to
### preweaning lethality terms

# import main files -------------------------------------------------------

### whole set of human protein coding genes with mouse othologues
### and mouse orthologues from hgnc
### there are some duplicates that we need to discard
### to keep only one2one relationships (strict criteria)

hgnc <- read_delim("https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt") %>%
  rename(gene_symbol = symbol, mgi_id = mgd_id) %>%
  select(hgnc_id, mgi_id) %>%
  separate_rows(mgi_id, sep = "\\|") %>%
  filter(!is.na(mgi_id))

hgnc_dups <- unique(hgnc$hgnc_id[duplicated(hgnc$hgnc_id)])
mgi_dups <- unique(hgnc$mgi_id[duplicated(hgnc$mgi_id)])

one2one <- hgnc %>%
  filter(!mgi_id %in% mgi_dups) %>%
  filter(!hgnc_id %in% hgnc_dups)

### whole set of mouse protein coding genes

mouse_genes <- read_delim("https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt") %>%
  select("MGI Accession ID", "Feature Type") %>%
  filter(`Feature Type` == "protein coding gene") %>%
  distinct()

mouse_proteincoding_genes <- unique(mouse_genes$`MGI Accession ID`)


### viability data impc
### lethal genes are easy to retrieve, there is a viability report
### three viability outcomes: lethal, subviable and viable

via_impc <- read_delim("http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz", delim = ",") %>%
  filter(is.na(Comment)) %>%
  rename(gene_symbol_mm  = "Gene Symbol",
         mgi_id = "Gene Accession Id",
         viability_impc = "Viability Phenotype HOMs/HEMIs") %>%
  select(mgi_id, viability_impc) %>%
  distinct()

### get lethal phenotypes from the MGI based on the set of lethal
### terms imported through the RData object ("lethal_terms")

via_mgi <- read_delim("https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt",
                      col_names = FALSE) %>%
  select(7, 5) %>%
  distinct() %>%
  rename(mgi_id = X7, mp_term = X5) %>%
  filter(mgi_id %in% mouse_proteincoding_genes) %>%
  mutate(mp_term_lethal = ifelse(mp_term %in% lethal_terms, "y", "n")) %>%
  select(mgi_id, mp_term_lethal) %>%
  distinct() %>%
  arrange(mgi_id, mp_term_lethal) %>%
  group_by(mgi_id) %>%
  summarise(mgi_lethal_term = paste0(unique(mp_term_lethal),
                                     collapse = "|")) %>%
  mutate(viability_mgi = ifelse(mgi_lethal_term == "n", "viable",
                                "lethal")) %>%
  select(mgi_id, viability_mgi) %>%
  distinct()

# create and export a single gene file ------------------------------------


### this file contains single gene knockout viability outcomes from
### the impc and mgi resources

viablity_impc_mgi <- one2one %>%
  left_join(via_impc) %>%
  left_join(via_mgi) %>%
  replace(is.na(.), "-")


write.table(viablity_impc_mgi, "./mouse_viability_data.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)

# modified files ----------------------------------------------------------

## recoding viability outcomes to generate a file with
## gene pairs and lethality annotations

viability_impc_mgi_recode <- viablity_impc_mgi %>%
  mutate(viability_impc = recode(viability_impc, lethal = "lethal",
                                 subviable = "lethal", viable = "viable")) %>%
  mutate(viability_both = paste0(viability_impc, viability_mgi)) %>%
  mutate(viability_both = ifelse(viability_both %in% c("viablelethal",
                                                       "lethalviable"),
                                 "discrepancy", viability_both)) %>%
  mutate(lethality_all = ifelse(grepl("viable", viability_both), "nonlethal",
                                ifelse(grepl("lethal", viability_both),
                                       "lethal", "-"))) %>%
  mutate(lethality_impc = ifelse(viability_impc == "viable", "nonlethal",
                                 viability_impc)) %>%
  rename(lethality_mgi = viability_mgi) %>%
  mutate(lethality_mgi = ifelse(lethality_mgi == "viable", "nonlethal",
                                lethality_mgi)) %>%
  select(hgnc_id, mgi_id, lethality_impc, lethality_mgi, lethality_all)

### pairs impc

impc_only <- viability_impc_mgi_recode %>%
  select(hgnc_id, lethality_impc) %>%
  filter(lethality_impc != "-")

impc_pairs <- t(combn(impc_only$hgnc_id, 2))

impc_pairs <- as.data.frame(impc_pairs)

names(impc_pairs) <- c("gene_a", "gene_b")

impc_pairs_lethality <- impc_pairs %>%
  left_join(impc_only, by = c("gene_a" = "hgnc_id"), multiple = "all") %>%
  rename(gene_a_impc = lethality_impc) %>%
  left_join(impc_only, by = c("gene_b" = "hgnc_id"), multiple = "all") %>%
  rename(gene_b_impc = lethality_impc)


### pairs mgi

mgi_only <- viability_impc_mgi_recode %>%
  select(hgnc_id, lethality_mgi) %>%
  filter(lethality_mgi != "-")

mgi_pairs <- t(combn(mgi_only$hgnc_id, 2))

mgi_pairs <- as.data.frame(mgi_pairs)

names(mgi_pairs) <- c("gene_a", "gene_b")

mgi_pairs_lethality <- mgi_pairs %>%
  left_join(mgi_only, by = c("gene_a" = "hgnc_id"), multiple = "all") %>%
  rename(gene_a_mgi = lethality_mgi) %>%
  left_join(mgi_only, by = c("gene_b" = "hgnc_id"), multiple = "all") %>%
  rename(gene_b_mgi = lethality_mgi)


### pairs combined

all <- viability_impc_mgi_recode %>%
  select(hgnc_id, lethality_all) %>%
  filter(lethality_all != "-")

all_pairs <- t(combn(all$hgnc_id, 2))

all_pairs <- as.data.frame(all_pairs)

names(all_pairs) <- c("gene_a", "gene_b")

all_pairs_lethality <- all_pairs %>%
  left_join(all, by = c("gene_a" = "hgnc_id"), multiple = "all") %>%
  rename(gene_a_all = lethality_all) %>%
  left_join(all, by = c("gene_b" = "hgnc_id"), multiple = "all") %>%
  rename(gene_b_all = lethality_all)


# export gene pair files ---------------------------------------------------

write.table(impc_pairs_lethality, "./mouse_impc_pairs_lethality.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(mgi_pairs_lethality, "./mouse_mgi_pairs_lethality.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(all_pairs_lethality, "./mouse_all_pairs_lethality.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
