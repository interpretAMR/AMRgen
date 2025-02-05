#' Generate a Stacked Bar Plot of MIC Values Colored by Gene Symbol for Each Antibiotic
#'
#' This function creates a stacked bar plot using `ggplot2`, where the x-axis 
#' represents MIC (Minimum Inhibitory Concentration) values (only using the 
#' numbers and ignoring the measurement signs so far), the y-axis represents 
#' their frequency, and the bars are colored by the gene symbols. 
#' The plot is faceted by antibiotic. 
#' @param pheno_data A data.frame containing columns "BioSample", "Antibiotic", 
#' "Measurement sign" (not used yet), "MIC (mg/L)", "Disk diffusion (mm)","Resistance phenotype".
#' BioSample is the ID used to match the entries in geno_data. 
#' MIC (Minimum Inhibitory Concentration) and Disk diffusion are numeric values, the other character strings.
#' @param geno_data A data.frame containing columns "Name", "Gene symbol", 
#' "Class" and "Subclass"  
#'
#' @param pathogen_mo A string result from as.mo() specifying what pathogen 
#' is used to select the breakpoints, e.g. 'B_STPHY_AURS' 
#' or group the data by different pathogens.
#' 
#' @param measure_method A string specifying the method used to measure antimicrobial resistance.
#' The default is 'MIC' for Minimum Inhibitory Concentration, the other option is 'DISK' for Disk Diffusion'.
#' 
#' @param breakpoint_source A string representing the source for breakpoint data. By default, it is set 
#' to the EUCAST standard from two years ago.
#' 
#' @param breakpoint_type A string specifying the type of breakpoint to use, default is 'human'.
#' 
#' @param color_by A string specifying the variable used to color the plot elements. Options include 
#' 'gene' (default) to color by gene, or 'class' to color by gene class. Any other strings also defaults to gene.
#' 
#' @param plot_title A string specifying the title of the plot. The default title is 
#' "Frequency distribution of MIC Value for different genotypes".
#' 
#' @param leg_pos A string specifying the position of the plot's legend. The default is 'bottom', 
#' but other options include 'top', 'left', or 'right'.
#' 
#' @param ncol_leg An integer specifying the number of columns for the legend. By default, it is 
#' calculated based on the number of unique gene symbols in `geno_data`, rounded to 10 items per column, but minimally 2.
#' This ensures the legend is well-organized when there are many gene symbols.
#' 
#' @return A ggplot2 object representing the stacked bar plot.
#' 
#' @details
#' - The function automatically groups the data by `Antibiotic`, `MIC (mg/L)`, and `Gene symbol` and calculates the frequency of occurrences.
#' - The plot uses the `viridis` color scale with the "turbo" palette by default. You can replace it with your custom palette if desired.
#'
#' @examples
#' AMRgen_barplot_mic_gen(pheno_data %>%
#'                          filter(Antibiotic %in% c('tetracycline','ceftaroline', 'levofloxacin')), 
#'                        geno_data, 
#'                        measure_method = 'MIC',
#'                        pathogen_mo = 'B_STPHY_AURS')
#' 
#' AMRgen_barplot_mic_gen(pheno_data %>%
#'                          filter(Antibiotic %in% c('erythromycin','penicillin', 'levofloxacin')), 
#'                        geno_data, 
#'                        color_by = 'class',
#'                        measure_method = 'DISK',
#'                        pathogen_mo = 'B_STPHY_AURS')
#'
#' @import ggplot2
#' @import dplyr
#' @import viridis
#' @import AMR
#' @export
#' 
AMRgen_barplot_mic_gen <- function(pheno_data, geno_data, 
                                   pathogen_mo,
                                   measure_method = 'MIC',
                                   breakpoint_source = paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2),
                                   breakpoint_type = 'human',
                                   color_by = 'gene', #alt.: 'class'
                                   
                                   plot_title = paste0("Frequency distribution of ",measure_method," Value for different genotypes"),
                                   leg_pos = 'bottom',
                                   ncol_leg = max(2, round(length(unique(geno_data$`Gene symbol`))/10))
                                   ){
  # todo: 
  # what about measurement signs? (ignoring for now)
  # add a lot more parameters (ggplot2)
  # shared frequency axis?
  # handle warnings

  # 1. get breakpoint data ----
  if (toupper(breakpoint_source) %in% clinical_breakpoints$guideline) {
    breakpoints <- toupper(breakpoint_source) # Replace with actual path
    message("Using ", breakpoint_source, " breakpoints.")
  } else {
    message("Invalid breakpoint source. Using default: " ,
         paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2))
    breakpoints <- paste('EUCAST', as.numeric(format(Sys.Date(), "%Y"))-2)
  }
  
  if (breakpoint_type %in% clinical_breakpoints$type) {
    breakpoint_type <- breakpoint_type # Replace with actual path
  } else {
    message("Invalid breakpoint type. Using default: human")
    breakpoint_type <- 'human'
  }
  
  if(!measure_method %in% c('MIC', 'DISK')){
    stop('Please choose one of the methods "MIC" or "DISK" for method.')
  }
  
  ab_bp <- data.frame('Antibiotic' = unique(pheno_data$Antibiotic)) %>%
    mutate(ab = as.character(as.ab(Antibiotic)))
           
  clinical_breakpoints_here <- clinical_breakpoints %>%
    filter(type == breakpoint_type,
           guideline == breakpoint_source,
           method == measure_method,
           mo == pathogen_mo,
           ab %in% ab_bp$ab)
  
  ab_bp <- ab_bp %>%
    left_join(clinical_breakpoints_here %>%
                select(ab, breakpoint_S, breakpoint_R, disk_dose),
              by = 'ab')
  
  # 2. prep data ----
  merged_data <- pheno_data %>%
    inner_join(geno_data, 
               by = c("BioSample" = "Name"),
               relationship = "many-to-many") %>% 
    # decide on variable to be plotted
    mutate(plot_var = ifelse(measure_method == 'MIC', 
                             `MIC (mg/L)`, 
                             `Disk diffusion (mm)`))

  # just for testing:
    # %>% slice_head(n=10000)
  
  # todo add disk breakpoint
  mic_summary <- merged_data %>%
    group_by(Antibiotic, plot_var, `Gene symbol`) %>%
    summarise(Frequency = n(),
              .groups = "drop") %>% # Count occurrences of MIC values per gene
    left_join(ab_bp, by=join_by('Antibiotic')) %>%
    left_join(geno_data %>% select(`Gene symbol`, Class), 
              by = "Gene symbol") %>%
    mutate(fill_col = case_when(
      color_by == 'class' ~ ifelse(is.na(Class), 'NA', Class), 
      TRUE ~`Gene symbol`))
  
  #3. plot ----
  ggplot(mic_summary, aes(x = plot_var, y = Frequency, fill = fill_col)) +
    geom_bar(stat = "identity", position = "stack",
             width = 0.6) + 
    facet_wrap(~ Antibiotic, scales = "free", ncol = NULL) +  # Create a panel for each antibiotic
    theme_minimal() +
    #exchange to our own color palette?
    scale_fill_viridis_d(option = "turbo") +  # Use "turbo" or other options like "viridis", "plasma"
    labs(
      title = plot_title,
      x = ifelse(measure_method == 'MIC', 
                 "MIC (mg/L)", 
                 "Disk diffusion (mm)"),
      y = "Frequency",
      fill = "Gene Symbol"
    ) +
    theme(
      legend.position = leg_pos,  # Move the legend to the top
      legend.title = element_blank(),
      legend.text = element_text(size = 8),  # Reduce legend text size
      legend.key.size = unit(0.4, "cm"),  # Reduce the size of legend keys
      axis.text.x = element_text(hjust = 1),
      strip.text = element_text(size = 10, face = "bold") # Customize panel titles
    ) +
    guides(fill = guide_legend(ncol = ncol_leg)) + # Arrange legend items into 4 columns
    geom_vline(data = ab_bp, aes(xintercept = `breakpoint_S`), 
               color = "black", linetype = "dashed") +  # Add breakpoint_S lines
    geom_vline(data = ab_bp, aes(xintercept = `breakpoint_R`), 
               color = "black", linetype = "dashed")   # Add breakpoint_R lines
}

# questions
# can we really expect the input data to look like this?
# should the user be able to paas all the parameters you theoretically could to ggplot2?


# example ----
# library(tidyverse)
# library(ggplot2)
# library(AMR)
# library(viridis) #might replace by own colors
# setwd('S:/OE/FG37/BotheC/AMRHackathon')
# pheno_data <- read_tsv(fs::path(getwd(), 'Data visualisation', 'astSA_matched.tsv'))
# geno_data <- read_tsv(fs::path(getwd(), 'Data visualisation', 'afpSA_matched.tsv'))
# AMRgen_barplot_mic_gen(pheno_data %>%
#                          filter(Antibiotic %in% c('tetracycline','ceftaroline', 'levofloxacin')), 
#                        geno_data, 
#                        measure_method = 'MIC',
#                        pathogen_mo = 'B_STPHY_AURS')
# 
# AMRgen_barplot_mic_gen(pheno_data %>%
#                          filter(Antibiotic %in% c('erythromycin','penicillin', 'levofloxacin')), 
#                        geno_data, 
#                        color_by = 'class',
#                        measure_method = 'DISK',
#                        pathogen_mo = 'B_STPHY_AURS')

