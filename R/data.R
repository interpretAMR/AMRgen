#' E. coli NCBI AST Example Data
#'
#' A subset of E. coli phenotype data from the NCBI AST browser.
#' @format ## `ecoli_ast_raw` A data frame with `r NROW(ecoli_ast_raw)` rows and `r NCOL(ecoli_ast_raw)` columns representing unprocessed data from the NCBI AST browser.
#' Columns include:
#' - `#BioSample`: Sample identifier.
#' - `Scientific name`: Species identifier.
#' - `Antibiotic`: Antibiotic name.
#' - `Testing standard`: Interpretation standard (EUCAST or CLSI).
#' - `Measurement sign`: Measurement sign (>, <, =, etc.) relating to MIC measurement.
#' - `MIC (mg/L)`: Minimum inhibitory concentration.
#' - `Disk diffusion (mm)`: Disk diffusion zone.
#' - `Resistance phenotype`: Resistance call (SIR) as submitted.
#' - ...: Additional metadata columns from the NCBI AST export.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast_raw"


#' E. coli NCBI AST Example Data, Re-interpreted with AMR Package
#'
#' A subset of E. coli phenotype data from the NCBI AST browser.
#' @format ## `ecoli_ast` A data frame with `r NROW(ecoli_ast)` rows and `r NCOL(ecoli_ast)` columns representing reinterpreted data from the NCBI AST browser.
#' Columns include:
#' - `id`: Sample identifier, imported from the `#BioSample` column in the raw input.
#' - `spp_pheno`: Species identifier, interpreted from `Scientific name` using `as.mo`, used to interpret `ecoff` and `pheno` columns.
#' - `drug_agent`: Antibiotic code, interpreted from `Antibiotic` using `as.ab`, used to interpret `ecoff` and `pheno` columns.
#' - `ecoff`: WT/NWT classification, interpreted using `as.sir`.
#' - `pheno`: S/I/R classification, interpreted using `as.sir`.
#' - `mic`: Minimum inhibitory concentration, formatted using `as.mic`, used to interpret `ecoff` and `pheno` columns.
#' - `disk`: Disk diffusion zone, formatted using `as.disk`, used to interpret `ecoff` and `pheno` columns.
#' - `guideline`: Interpretation guidelines used to interpret `ecoff` and `pheno` columns.
#' - `Scientific name`: Species identifier, from input file.
#' - `Antibiotic`: Antibiotic name, from input file.
#' - `Testing standard`: Interpretation standard (EUCAST or CLSI), from input file.
#' - `Measurement sign`: Measurement sign (>, <, =, etc.) relating to MIC measurement, from input file.
#' - `MIC (mg/L)`: Minimum inhibitory concentration, from input file.
#' - `Disk diffusion (mm)`: Disk diffusion zone, from input file.
#' - `Resistance phenotype`: Resistance call (SIR), from input file.
#' - ...: Additional metadata columns from the NCBI AST export.
#' @source <https://www.ncbi.nlm.nih.gov/pathogens/ast>
"ecoli_ast"


#' E. coli Genotype Example Data
#'
#' Genotypes called using AMRfinderPlus (v3.12.8, DB 2024-01-31.1), sourced from the AllTheBacteria project.
#' @format ## `ecoli_geno_raw` A data frame with `r NROW(ecoli_geno_raw)` rows and `r NCOL(ecoli_geno_raw)` columns representing genotyping results from AMRfinderPlus.
#' Columns include:
#' - `Name`: Sample identifier.
#' - `Gene symbol`: Gene symbol in NCBI RefGene.
#' - `Hierarchy node`: Node in NCBI hierarchy.
#' - `Class`, `Subclass`: Drug class(es) associated with the marker (from NCBI RefGene).
#' - `% Coverage of reference sequence`, `% Identity to reference sequence`, `Accession of closest sequence`: Sequence match information.
#' - ...: Additional metadata columns from the AMRfinderPlus output.
#' @source <https://github.com/ncbi/amr/wiki/Interpreting-results>
"ecoli_geno_raw"
