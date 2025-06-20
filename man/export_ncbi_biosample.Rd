% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export_ncbi.R
\name{export_ncbi_biosample}
\alias{export_ncbi_biosample}
\title{Import/Export BioSample Antibiograms}
\usage{
export_ncbi_biosample(data, file, overwrite = FALSE)
}
\arguments{
\item{data}{Data set containing SIR results.}

\item{file}{File path to which the data will be exported.}

\item{overwrite}{Logical indicating whether to overwrite an existing file.}
}
\value{
A tab-delimited UTF-8 text file is written to disk in the format required by NCBI BioSample Antibiograms.
}
\description{
Output phenotype data to \href{https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/}{NCBI BioSample Antibiograms}.
}
