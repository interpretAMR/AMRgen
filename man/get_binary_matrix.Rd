% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_binary_matrix.R
\name{get_binary_matrix}
\alias{get_binary_matrix}
\title{Get Binary Matrix of Genotype and Phenotype Data}
\usage{
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic,
  drug_class_list,
  keep_SIR = TRUE,
  keep_assay_values = FALSE,
  keep_assay_values_from = c("mic", "disk"),
  geno_sample_col = NULL,
  pheno_sample_col = NULL,
  sir_col = "pheno",
  ecoff_col = "ecoff",
  marker_col = "marker",
  most_resistant = TRUE
)
}
\arguments{
\item{geno_table}{A data frame containing genotype data, including at least one column labeled \code{drug_class} for drug class information and one column for sample identifiers (specified via \code{geno_sample_col}, otherwise it is assumed the first column contains identifiers).}

\item{pheno_table}{A data frame containing phenotype data, which must include a column \code{drug_agent} (with the antibiotic information), a column with the resistance interpretation (S/I/R, specified via \code{sir_col}), and optionally a column with the ECOFF interpretation (WT/NWT, specified via \code{ecoff_col}).}

\item{antibiotic}{A character string specifying the antibiotic of interest to filter phenotype data. The value must match one of the entries in the \code{drug_agent} column of \code{pheno_table}.}

\item{drug_class_list}{A character vector of drug classes to filter genotype data for markers related to the specified antibiotic. Markers in \code{geno_table} will be filtered based on whether their \code{drug_class} matches any value in this list.}

\item{keep_SIR}{A logical indicating whether to retain the full S/I/R phenotype column in the output. Defaults to \code{TRUE}.}

\item{keep_assay_values}{A logical indicating whether to include columns with the raw phenotype assay data. Assumes there are columns labelled "mic" and "disk"; these will be added to the output table if present. Defaults to \code{FALSE}.}

\item{keep_assay_values_from}{A character vector specifying which assay values (e.g., \code{"mic"}, \code{"disk"}) to retain if \code{keep_assay_values} is \code{TRUE}. Defaults to \code{c("mic", "disk")}.}

\item{geno_sample_col}{A character string (optional) specifying the column name in \code{geno_table} containing sample identifiers. Defaults to \code{NULL}, in which case it is assumed the first column contains identifiers.}

\item{pheno_sample_col}{A character string (optional) specifying the column name in \code{pheno_table} containing sample identifiers. Defaults to \code{NULL}, in which case it is assumed the first column contains identifiers.}

\item{sir_col}{A character string specifying the column name in \code{pheno_table} that contains the resistance interpretation (SIR) data. The values should be interpretable as "R" (resistant), "I" (intermediate), or "S" (susceptible).}

\item{ecoff_col}{A character string specifying the column name in \code{pheno_table} that contains the ECOFF interpretation of phenotype. The values should be interpretable as "WT" (wildtype) or "NWT" (nonwildtype).}

\item{marker_col}{A character string specifying the column name in \code{geno_table} containing the marker identifiers. Defaults to \code{"marker"}.}

\item{most_resistant}{A logical indicating whether, when multiple phenotype entries are present for the same sample and drug, the most resistant should be kept (otherwise the least resistant is kept). Default is \code{TRUE}.}
}
\value{
A data frame where each row represents a sample, and each column represents a genetic marker related to the specified antibiotic's drug class. The binary values in the matrix indicate the presence (\code{1}) or absence (\code{0}) of each marker for each sample, along with resistance status columns for the specified antibiotic: \code{R} for resistant (defined from \code{sir_col}, 1=R, 0=I/S) and \code{NWT} for nonwildtype (defined by \code{ecoff_col} if provided: 1=NWT, 0=WT; otherwise defined from \code{sir_col}: 1=I/R, 0=S).
}
\description{
This function generates a binary matrix representing the resistance (R vs S/I) and nonwildtype (R/I vs S) status for a given antibiotic, and presence or absence of genetic markers related to one or more specified drug classes. It takes as input separate tables for genotype and phenotype data, matches these according to a common identifier (either specified by column names or assuming the first column contains the ID), and filters the data according to the specified antibiotic and drug class criteria before creating a binary matrix. Suitable input files can be generated using \code{import_ncbi_ast} to import phenotype data from NCBI, and \code{parse_amrfp} to import genotype data from AMRfinderPlus.
}
\details{
This function performs several steps:
\itemize{
\item Verifies that the \code{pheno_table} contains a \code{drug_agent} column and converts it to class \code{ab} if necessary.
\item Filters the \code{pheno_table} to retain data related to the specified antibiotic.
\item Checks that the \code{geno_table} contains markers associated with the specified drug class(es).
\item Matches sample identifiers between \code{geno_table} and \code{pheno_table}.
\item Extracts and transforms the phenotype data into a binary format indicating resistance and NWT status.
\item Constructs a binary matrix for genotype data, with each column representing a genetic marker.
\item Returns a single matrix with sample identifiers plus binary variables for each phenotype and genetic marker.
}
}
\examples{
\dontrun{
geno_table <- parse_amrfp("testdata/Ecoli_AMRfinderplus_n50.tsv", "Name")
pheno_table <- import_ncbi_ast("testdata/Ecoli_AST_NCBI_n50.tsv")
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype"
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE,
  keep_assay_values_from = "mic"
)
get_binary_matrix(
  geno_table,
  pheno_table,
  antibiotic = "Ciprofloxacin",
  drug_class_list = c("Quinolones"),
  sir_col = "Resistance phenotype",
  keep_assay_values = TRUE,
  keep_assay_values_from = "MIC (mg/L)"
)
}
}
\seealso{
\code{compare_geno_pheno_id}, \code{as.ab}, \code{as.sir}, \code{ab_name}
}
