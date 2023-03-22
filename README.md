### Data sources:
-  Mouse viability data IMPC: IMPC DR 18.0 http://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/release-18.0/results/
-  Mouse viability data MGI: https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt"
-  Set of lethal terms (Mammalian Phneotype Ontology (MP) from Dickinson et al. 2016 https://pubmed.ncbi.nlm.nih.gov/27626380/
-  Human - mouse orthologues: "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt
-  Date of access: 22/03/2023

### Notes:
-  HGNC ids are used as gene stable identifiers
-  Annotations for human / mouse protein coding genes
-  Only one2one orthologues from the HGNC file were kept
-  IMPC mouse viability: standardised screens, primary viability assessment: 3 outcomes (lethal, subviable, viable). These categories are simplified here as 'lethal' (lethal % subviable) and 'nonlethal' (viable) - lethality from the embryonic period to preweaning age).
-  MGI mouse viability: MGI data contains data collated from the literature, from heterogeneous lines, it includes data from the IMPC resource. The genes labelled as 'lethal' are those genes with a mouse knockout with an associated mammalian phenotype ontology term from a curated set of terms. The genes labelled as 'viable' include any gene where none of the mouse models included in the resource are associated to any of those terms
-  All mouse viability: IMPC and MGI data combined. Any discrepancies in lethality between the 2 resources are excluded.
-  **Gene pairs are unique pairs, e.g. gene_a - gene_b won't show as gene_b - gene_a**
-  **Data on mouse lethality come from single gene knockouts, these are not synthetic lethality pairs**.