#' PlotDEGenes
#'
#' Using a bubbleplot, plot genes that are differentially expressed between
#' cells of a comparison group within a grouping (i.e. DEGs between different
#' disease groups within the same cell type).  The results of all intragroup
#' comparisons are compiled together and plotted together on one bubbleplot.
#'
#' @param seuratObj Processed Seurat scRNAseq object
#' @param combo_list Differentially expressed gene list for combination of
#'   contrast factors from FindGroupDEGs
#' @param group_by Group variable within which to examine the cells.
#'   Default: celltype
#' @param compare_by Group variable to use when comparing the cells.
#'   Default: class.
#' @param ... Extra arguments to pass to seuratBubblePlot::bubbleplot()
#'
#' @import magrittr
#' @importFrom future plan multiprocess
#' @importFrom glue glue
#' @importFrom Seurat SubsetData SetAllIdent
#' @importFrom future.apply future_lapply
#' @importFrom seuratBubblePlot bubbleplot
#' @importFrom stringr str_wrap
#' @importFrom plyr mapvalues
#'
#' @return ggplot2 object
#' @export PlotDEGenes
#'
#' @examples
PlotDEGenes <- function(seuratObj,
                        combo_list,
                        group_by = "celltype",
                        compare_by = "class",
                        gene_annotation_df = NULL,
                        split_if = NULL,
                        ...) {
    seuratObj <- SetAllIdent(seuratObj, group_by)
  # for each cell type
  opl <- lapply(
    names(combo_list),
    function(x) {
      if (length(combo_list[x][[1]]) > 0) {
        print(glue("Now plotting: {x}"))
        celltypeObj <- SubsetData(seuratObj,
          ident.use = x,
          subset.raw = TRUE
        )
        celltypeObj <- SetAllIdent(
          celltypeObj,
          compare_by
        )
        common_genes <-
          unique(unlist(lapply(
            combo_list[[x]],
            function(y) {
              y[, "gene"]
            }
          )))

        if (length(common_genes) > 0) {
          if (!is.null(gene_annotation_df)) {
            common_genes %<>% as.data.frame()
            colnames(common_genes) <- "genes"
            common_genes$annotations <- mapvalues(
              x = common_genes$genes,
              from = gene_annotation_df$genes,
              to = gene_annotation_df$annotations
            )
            annotated_gene_list <- TRUE
          } else {
            annotated_gene_list <- FALSE
          }


          bp <- bubbleplot(
                            celltypeObj,
                            genes_plot = common_genes,
                            x_axis_title = str_wrap(x, 20),
                            y_axis_title = NULL,
                            do_return = TRUE,
                            pct_legend_title = "Percent\ngroup\nexpressing",
                            scale_legend_title = "Average\nscaled\nexpression",
                            annotated_gene_list = annotated_gene_list,
                            ...
                          )
          if(!is.null(split_if) & (length(common_genes) > split_if)){
            d <- length(unique(bp$data$genes_plot))
            split_groups <- rep(1:ceiling(d/20), each = split_length)[1:length(unique(bp$data$genes_plot))]
            split_table <-  data.frame('genes' = unique(bp$data$genes_plot),
                                       'split_group' = split_groups)

            alpha$data$split_groups <- mapvalues(x = bp$data$genes_plot,
                                                 from = split_table$genes,
                                                 to = split_table$split_group)
            bp + facet_grid(.~split_group, scales = "free") + theme(strip.text = element_none())
          } else {
            bp
          }
        } else {
          print(glue("No genes to plot for {x} were found."))
        }
      }
    }
  )
  return(opl)
}
