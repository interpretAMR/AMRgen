### Function to generate "amr_upset"

#' amr_upset: Generate a series of plots for AMR gene and combination analysis
#'
#' This function generates a set of visualizations to analyze AMR gene combinations, MIC values, and gene prevalence
#' from a given binary matrix. It creates several plots, including MIC distributions, a bar plot for
#' the percentage of strains per combination, a dot plot for gene combinations in strains, and a plot for gene prevalence.
#'
#' @param binary_matrix A data frame containing the original binary matrix output from the `get_binary_matrix` function,
#'        with columns representing genes, resistance, MIC values, and metadata such as microorganism and antibiotic
#'        information. This needs to be updated / standardised in future versions of AMRgen.
#' @param min_set_size An integer specifying the minimum size for a gene set to be included in the analysis and plots.
#'        Default is 2. Only genes with at least this number of occurrences are included in the plots.
#' @param order A character string indicating the order of the combinations on the x-axis. Options are:
#'        - "": Default (decreasing frequency of combinations),
#'        - "genes": Order by the number of genes in each combination,
#'        - "mic": Order by the median MIC of each combination. Default is decreasing frequency.
#' @param plot_set_size Logical indicating whether to include a bar plot showing the set size (i.e. 
#'        number of times each combination of markers is observed). Default is FALSE.
#' @param plot_category Logical indicating whether to include a stacked bar plot showing, for each marker combination,
#'         the proportion of samples with each phenotype classification (specified by the 'pheno' column in the input file). 
#'         Default is TRUE.
#'
#' @return A list containing the following elements:
#'   \describe{
#'     \item{plot}{A grid of plots displaying: (i) grid showing the marker combinations observed, MIC distribution per marker combination, frequency per marker and (optionally) phenotype classification and/or number of samples for each marker combination.}
#'     \item{mic_summary}{Summary of each marker combination observed, including median MIC (and interquartile range) and positive predictive value for resistance (R).}
#'   }
#'
#' @details This function processes the provided binary matrix (`binary_matrix`), which is expected to contain data on gene 
#'          combinations found in the strain, MIC for that strain (e.g., resistance or susceptibility) and 
#'          corresponding MIC values for different genes.
#'          The function also includes an analysis of gene prevalence and an ordering option for visualizing combinations
#'          in different ways.
#'
#' @examples
#' \dontrun{
#' # Example usage
#'
#' ecoli_geno <- import_amrfp(ecoli_geno_raw, "Name")
#' 
#' binary_matrix<- get_binary_matrix(geno_table=ecoli_geno, 
#'               pheno_table=ecoli_ast, 
#'               antibiotic="Ciprofloxacin", 
#'               drug_class_list=c("Quinolones"), 
#'               sir_col="pheno", 
#'               keep_assay_values=TRUE, 
#'               keep_assay_values_from = "mic"
#'            )
#' 
#' amr_upset(binary_matrix, min_set_size = 3, order = "mic")
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom forcats fct_rev
#' @importFrom AMR as.mic
#' @import patchwork
#' 
#' @export
amr_upset <- function(binary_matrix, min_set_size = 2, order = "", 
                      plot_set_size=FALSE, plot_category=TRUE) {
  ## Inputs
  # takes in binary_matrix = output from get_binary_matrix function
  # takes in order = single value. Default is decreasing frequency. 
  #         "genes" = # genes. "mic" = median mic 
  # default min set size is 2 (greater than or equal to this)
  
  # tidy up binary_matrix
  col <- colnames(binary_matrix) # get column names 
  
  # extract only the gene column names - need to exclude mic, disk, R, NWT (standard col names)
  # and the id column which will be the first col, doesn't matter what it's called
  # remaining columns will be the genes
  cols_to_remove <- c("mic", "disk", "R", "NWT", "pheno")
  genes <- col[-1]

  # gene names
  genes <- setdiff(genes, cols_to_remove)

  # Add in a combination column and filter to samples with MIC data only
  binary_matrix_wide <- binary_matrix %>% 
    filter(!is.na(mic)) %>% # make this optionally disk also
    unite("combination_id", genes[1]:genes[length(genes)], remove = FALSE)  # add in combinations 
  
  # Make matrix longer 
  binary_matrix <- binary_matrix_wide %>% pivot_longer(cols = genes[1]:genes[length(genes)], names_to = "genes")
  
  ##### Data wrangling for plots ###
  ### Bar plot - X axis = combination. Y axis = number of strains ###
  # This first to filter on combinations with enough data 
  bar_plot <- binary_matrix_wide %>% group_by(combination_id) %>%
    summarise(n = n()) %>%
    mutate(perc = 100 * n / sum(n)) %>% # count number with each combination
    filter(n > 1) # count filter
  
  # which have enough strains / data
  comb_enough_strains <- bar_plot %>% pull(combination_id)

  ### MIC plot - dot plot. X axis = combination. Y axis = MIC #####
  mic_plot <- binary_matrix %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(combination_id, mic, R) %>%
    summarise(n = n()) # count how many at each MIC, keep resistant for colour
  
  ### Gene prevalence plot 
  gene.prev <- binary_matrix %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(genes) %>%
    summarise(gene.prev = sum(value))

  ################### ORDER y axis ########################
  ### Set order of genes for y axis in combination dot plot 
  ## And filter out on set size mininum
  gene.order.desc <- gene.prev %>%
    arrange(desc(gene.prev)) %>%
    filter(gene.prev >= min_set_size) %>%
    pull(genes)

  # For gene prev plot
  gene.prev <- gene.prev %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc))
  
  # For combination dot plot 
  binary_matrix <- binary_matrix %>%
    filter(combination_id %in% comb_enough_strains) %>% 
    filter(genes %in% gene.order.desc) %>% 
    mutate(genes = factor(genes, levels = gene.order.desc))
  

  
  ############# Which have lines between in dot plot? 
  ### Point plot - X axis = combination. Y axis = genes. Lines joining genes in same strain ### 
  binary_matrix$point_size = 2 * binary_matrix$value # want dot size to be larger than 1 => can make 2/3/4 etc
  # Only plot a line between points if more than one gene in a strain 
  # get how many in each strain 
  multi_genes_combination_id_all <- binary_matrix %>%
    group_by(combination_id) %>%
    filter(value == 1) %>%
    mutate(u = length(unique(genes))) %>%
    filter(genes %in% gene.order.desc) %>%
    mutate(genes = factor(genes, levels = gene.order.desc))
  # get only those with > 1
  multi_genes_combination_ids <- multi_genes_combination_id_all %>%
    filter(u > 1) %>%
    mutate(
      min = first(genes),
      max = last(genes)
    ) %>%
    ungroup()

  ################### ORDER x axis ########################


  ### Set order of combination_id <- x axis
  # Default = decreasing frequency
  ordered_comb_order <- bar_plot %>%
    arrange(desc(perc)) %>%
    pull(combination_id)
  mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
  bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
  binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  
  # Do by # genes in combination (only want each id once)
  if (order == "genes") {
    ordered_comb_order <- multi_genes_combination_id_all %>%
      arrange(u) %>%
      filter(row_number() == 1) %>%
      pull(combination_id)
    mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
    bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
    binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  }
  # Do by # median mic in combination (only want each id once)
  if(order == "mic"){
    mic_medians <- binary_matrix_wide %>% group_by(combination_id) %>% summarise(median = median(mic))
    ordered_comb_order <- mic_medians %>% arrange(median) %>% pull(combination_id)
    mic_plot$combination_id <- factor(mic_plot$combination_id, levels = ordered_comb_order)
    bar_plot$combination_id <- factor(bar_plot$combination_id, levels = ordered_comb_order)
    binary_matrix$combination_id <- factor(binary_matrix$combination_id, levels = ordered_comb_order)
  }
  
  ##### Plots ###
  ### AMR package colours
  colours_SIR <- c("#3CAEA3", "#F6D55C", "#ED553B")
  # currently only 0/1
  # names(colours_SIR) <- c("S", "I", "R")

  ### MIC plot
  g1 <- ggplot(data = mic_plot, aes(combination_id, mic)) +
    geom_point(aes(size = n, colour = as.factor(R)), show.legend = TRUE) +
    # geom_hline(data = cut_dat, aes(yintercept = AMR::as.mic(breakpoint_S)), colour = colours_SIR["S"]) +
    # geom_hline(data = cut_dat, aes(yintercept = AMR::as.mic(breakpoint_R)), colour = colours_SIR["R"]) +
    theme_bw() +
    scale_y_mic() +
    ylab("MIC (mg/L)") +
    scale_x_discrete("group") +
    scale_size_continuous("Number of\nisolates") +
    scale_color_manual("Resistance\nclass", values = colours_SIR) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  ### Bar plot  
  g2 <- ggplot(bar_plot, aes(x=combination_id, y = perc)) + 
    geom_bar(stat = "identity") + 
    theme_bw() +
    #geom_text(aes(label = n), nudge_y = -.5) + 
    scale_y_reverse("Percentage") +
    #scale_y_continuous("Percentage") + 
    scale_x_discrete("group") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
  
  category_plot <- binary_matrix %>% select(id, pheno, combination_id) %>% distinct() %>%
    ggplot(aes(x=combination_id, fill=pheno)) +
    geom_bar(stat = "count", position = "fill") +
    #scale_fill_manual(values = plot_cols) +
    #geom_text(aes(label = ..count..), stat = "count", position = position_fill(vjust = .5), size = 3) +
    theme_light() +
    labs(x = "", y = "Category", fill = "Category") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
    
  ### Dot plot of combinations 
  g3 <- binary_matrix %>% 
    mutate(binary_comb=if_else(value>0, 1, 0)) %>%
    ggplot(aes(combination_id, fct_rev(genes))) + 
    geom_point(aes(size = binary_comb), show.legend = FALSE) + 
    theme_bw() + 
    scale_size_continuous(range = c(-1,2)) + 
    scale_x_discrete("group") + 
    scale_y_discrete("Marker") + 
    geom_segment(data = multi_genes_combination_ids,
                 aes(x = combination_id, xend = combination_id, 
                     y = min, yend = max, group = combination_id),
                 color = "black") + # add lines => dashed ok? 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

  ### Plot gene prev / set size
  g4 <- ggplot(gene.prev, aes(x = fct_rev(genes), y = gene.prev)) +
    geom_col() +
    theme_bw() +
    coord_flip() +
    scale_y_reverse("Set size") +
    theme(axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank())
  
  # assemble plot
  final_plot <- patchwork::plot_spacer() + g1 + patchwork::plot_layout(ncol = 2, widths = c(1, 4), guides="collect")
  if (plot_category) { final_plot <- final_plot + patchwork::plot_spacer() + category_plot }
  final_plot <- final_plot + g4 + g3
  if (plot_set_size) { final_plot <- final_plot + patchwork::plot_spacer() + g2 }

  # set relative plotting heights
  if (plot_category & plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,1,2,1))}
  if (plot_category & !plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,1,2))}
  if (!plot_category & plot_set_size) {final_plot <- final_plot + patchwork::plot_layout(heights=c(2,2,1))}
  
  print(final_plot)
  
  # summary table
  mic_summary <- binary_matrix_wide %>% 
    filter(combination_id %in% comb_enough_strains) %>% 
    group_by(combination_id) %>%
    summarise(median = median(mic), 
              lower=stats::quantile(mic,0.25),
              upper=stats::quantile(mic,0.75),
              ppv=mean(R, na.rm=T),
              R=sum(R, na.rm=T),
              n=n())
  
  # get names for mic_summary
  combination_names <- binary_matrix_wide %>% 
    select(combination_id, any_of(genes)) %>% 
    distinct() %>% 
    filter(combination_id %in% comb_enough_strains)
  
  combination_names <- combination_names %>% 
    mutate(marker_list = apply(., 1, function(row) {
      paste(names(combination_names)[-1][row[-1] == 1], collapse = ", ")
    }), .after=combination_id) %>%
    mutate(marker_count = rowSums(. == 1)) %>%
    select(combination_id, marker_list, marker_count)
  
  mic_summary <- mic_summary %>% 
    left_join(combination_names) %>%
    select(-combination_id) %>%
    mutate(marker_list=if_else(is.na(marker_list), "-", marker_list)) %>%
    mutate(marker_count=if_else(is.na(marker_count), 0, marker_count)) %>% 
    relocate(marker_list, .before=median) %>% 
    relocate(marker_count, .before=median)
  
  return(list(plot=final_plot, mic_summary=mic_summary))
}
