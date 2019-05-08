#' @title PlotDEGenes
#'
#' @description Using a bubbleplot, plot genes that are differentially expressed between
#' cells of a comparison group within a grouping (i.e. DEGs between different
#' disease groups within the same cell type).  The results of all intragroup
#' comparisons are compiled together and plotted together on one bubbleplot.
#'
#' @param object Processed Seurat scRNAseq object
#' @param combo_list Differentially expressed gene list for combination of
#' contrast factors from FindGroupDEGs
#' @param group_by Group variable within which to examine the cells.
#' Default: celltype
#' @param compare_by Group variable to use when comparing the cells.
#' Default: class.
#' @param gene_annotation_df Dataframe containing
#' @param split_if If the number of features is above this amount, facet the plot. Default: NULL
#' @param ... Extra arguments to pass to SeuratBubblePlot::bubbleplot()
#'
#' @importFrom glue glue
#' @importFrom Seurat SubsetData Idents<-
#' @importFrom future.apply future_lapply
#' @importFrom SeuratBubblePlot bubbleplot
#' @importFrom stringr str_wrap
#' @importFrom dplyr recode pull
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom ggplot2 facet_grid theme element_blank
#'
#' @return ggplot2 object
#' @export PlotDEGenes
#'
#' @examples
PlotDEGenes <- function(object,
                        combo_list,
                        group_by = "celltype",
                        compare_by = "class",
                        gene_annotation_df = NULL,
                        split_if = NULL,
                        ...) {
  Suggested.Symbol <- NULL

  Idents(object) <- object[[group_by]]
  # for each cell type
  opl <- lapply(
    names(combo_list),
    function(x) {
      if (length(combo_list[x][[1]]) > 0) {
        message(glue("Now plotting: {x}"))
        celltypeObj <- SubsetData(object,
                                  ident.use = x)
        Idents(celltypeObj) <- celltypeObj[[compare_by]]
        common_genes <- lapply(combo_list[[x]], function(y){ y[,"gene"]}) %>%
          unlist() %>%
          unique()

        if (length(common_genes) > 0) {
          if (!is.null(gene_annotation_df)) {
            common_genes %<>% as.data.frame()
            colnames(common_genes) <- "genes"

            annotations <- gene_annotation_df[["annotations"]]
            names(annotations) <- gene_annotation_df[["genes"]]

            common_genes$annotations <- recode(common_genes[["genes"]],
                                               !!!annotations)
            annotated_gene_list <- TRUE
          } else {
            annotated_gene_list <- FALSE
          }

          common_genes %<>%
            checkGeneSymbols() %>%
            filter(!is.na(Suggested.Symbol)) %>%
            pull(Suggested.Symbol)

          bp <- bubbleplot(celltypeObj,
                           genes_plot = common_genes,
                           x_axis_title = str_wrap(x, 20),
                           y_axis_title = NULL,
                           do_return = TRUE,
                           pct_legend_title = "Percent\ngroup\nexpressing",
                           scale_legend_title = "Average\nscaled\nexpression",
                           annotated_gene_list = annotated_gene_list,
                           translate_gene_names = TRUE,
                            ...)
          if (!is.null(split_if)){
            if (length(common_genes) > split_if){
              d <- bp$data$genes_plot %>% unique() %>% length()
              split_groups <- rep(1:ceiling(d/20), each = split_if)[1:d]
              split_table <-  data.frame('genes' = unique(bp$data$genes_plot),
                                         'split_group' = split_groups)

              new_split_group <- split_table[["split_group"]]
              names(new_split_group) <- split_table[["genes"]]
              alpha$data$split_groups <- recode(bp$data$genes_plot,
                                                   !!!new_split_group)
              bp + facet_grid(.~split_group, scales = "free") + theme(strip.text = element_blank())
            }
          } else {
            bp
          }
        } else {
          message(glue("No genes to plot for {x} were found."))
        }
      }
    }
  )
  return(opl)
}
