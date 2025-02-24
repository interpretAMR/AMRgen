# # hAMRonize function & test code
# # Zoe A. Dyson (zoe.dyson@lshtm.ac.uk)
# # Last updated 14/01/2025
#
# # List required packages
# packages <- c("reticulate", "tidyverse")
#
# # Install missing packages
# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }
#
# # Load packages
# invisible(lapply(packages, library, character.only = TRUE))
#
# # create a new environment for python packages required for hAMRonization
# conda_create("amrgen-reticulate")
#
# # install hAMRonization python package
# conda_install(
#   envname = "amrgen-reticulate",
#   packages = "hAMRonization",
#   forge = T,
#   pip = T
# )

# # select python environment
# use_condaenv("amrgen-reticulate")

#' @title hamronize_data
#'
#' @param user_software_name the analysis software used to screen genome data for AMR determinants (must be amrfinderplus, rgi, or resfinder)
#' @param user_software_version the version of the analysis software used to screen genome data
#' @param user_database_version the version of the database used
#' @param user_input_filename the name of the genotypic AMR data file
#'
#' @return A data frame containing 'harmonized' AMR genotype data
#'
#' @examples
#' \dontrun{
#' harmonize_data(
#'   user_software_name = "amrfinderplus",
#'   user_software_version = "3.12.8",
#'   user_input_filename = "ATB_Achromobacter_AFP.tsv",
#'   user_database_version = "2024-01-31.1"
#' )
#' }
harmonize_data <- function(user_software_name,
                           user_software_version,
                           user_database_version,
                           user_input_filename) {
  # send user input data to python
  r_to_py(user_input_filename, convert = T)
  r_to_py(user_software_name, convert = T)
  r_to_py(user_software_version, convert = T)
  r_to_py(user_database_version, convert = T)

  # import python libraries
  py_run_string("import hAMRonization")
  py_run_string("import pandas") # might not need to call explicitly as hamronization dependency

  # create empty list to store hamronized results
  py_run_string("hamronized_output = []")

  # make python variables from R variables (required for table filling)
  py_run_string("filename = r.user_input_filename")
  py_run_string("software_name = r.user_software_name")
  py_run_string("software_version = r.user_software_version")
  py_run_string("database_version = r.user_database_version")

  # run hamronization
  py_run_string("metadata = {'analysis_software_version': software_version, 'reference_database_version': database_version, 'input_file_name': filename}")
  py_run_string("for result in hAMRonization.parse(filename, metadata, software_name): hamronized_output.append(result)")

  # convert output to data frame
  py_run_string("hamronized_output_df = pandas.DataFrame(hamronized_output)")

  # convert output to AMRgen genotype data frame format
  hamronized_data <- py$hamronized_output_df %>%
    mutate(Sample_ID = input_sequence_id) %>%
    mutate(marker = as.gene(gene_symbol)) %>%
    mutate(drug_agent = antimicrobial_agent) %>%
    select(Sample_ID, marker, drug_class, drug_agent)
  
  # Separate drug classes for rgi data
  if (user_software_name=="rgi"){
    hamronized_data <- hamronized_data %>%
      separate_longer_delim(., drug_class, delim=";") %>%
      mutate(drug_class = str_trim(drug_class, side = "both"))
  }
  
  # return harmonized data frame
  return(hamronized_data)
}
