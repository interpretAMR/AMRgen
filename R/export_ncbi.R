#' Import/Export BioSample Antibiograms
#'
#' Output phenotype data to [NCBI BioSample Antibiograms](https://www.ncbi.nlm.nih.gov/biosample/docs/antibiogram/).
#' @param data Data set containing SIR results.
#' @param file File path to which the data will be exported.
#' @param overwrite Logical indicating whether to overwrite an existing file.
#' @importFrom AMR is.sir
#' @importFrom dplyr transmute where
#' @importFrom tidyr pivot_longer
#' @return A tab-delimited UTF-8 text file is written to disk in the format required by NCBI BioSample Antibiograms.
#' @export
export_ncbi_biosample <- function(data, file, overwrite = FALSE) {
  if (file.exists(file) && !overwrite) {
    stop("The file ", file, " already exists and `overwrite`` is set to `FALSE`")
  }
  if (!grepl("[.](txt|tsv)$", file, ignore.case = TRUE)) {
    stop("`file` must have the file extension 'txt' or 'tsv'")
  }
  if (any(is.sir(data))) {
    data <- data %>%
      pivot_longer(where(AMR::is.sir),
        names_to = "antibiotic",
        values_to = "resistance_phenotype"
      )
  }
  out <- data %>%
    transmute(
      `sample_name/biosample_accession` = NA,
      antibiotic = NA,
      resistance_phenotype = NA,
      measurement_sign = NA,
      measurement = NA,
      measurement_units = NA,
      laboratory_typing_method = NA,
      laboratory_typing_platform = NA,
      vendor = NA,
      laboratory_typing_method_version_or_reagent = NA,
      testing_standard = NA
    )
  utils::write.table(
    x = out,
    file = file,
    append = FALSE,
    quote = TRUE,
    sep = "\t",
    na = "",
    row.names = FALSE,
    col.names = TRUE,
    fileEncoding = "UTF-8"
  )
}
