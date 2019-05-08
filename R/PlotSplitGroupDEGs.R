#' PlotSplitGroupDEGs
#'
#' For a given Seurat object and a list of intragroup DEGs, for each grouping
#' plot the genes that were differentially expressed between each comparison_group.
#' The results for each comparison within the group will be plotted in its own
#' bubbleplot, and all comparisons for that group are combined into one compound plot.
#'
#' @param object Processed Seurat scRNAseq object
#' @param combo_list List of list of differentially expressed genes from DEcombinations
#' @param group_by Primary grouping variable (from ident slot or meta.data column)
#' @param compare_by Intragroup comparison variable
#' @param ... Additional parameters to pass to bubbleplot
#'
#' @return A list of gtables
#' @export PlotSplitGroupDEGs
#'
#' @importFrom glue glue
#' @importFrom Seurat Idents<- SubsetData
#' @importFrom SeuratBubblePlot bubbleplot
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggplot2 theme
#' @importFrom cowplot get_legend
#'
#' @examples
PlotSplitGroupDEGs <- function(object,
                               combo_list,
                               group_by = "celltype",
                               compare_by = "class",
                               ...) {
  Idents(object) <- object[[group_by]]
  # for each cell type
  opl <- lapply(
    names(combo_list),
    function(x) {
      if (length(combo_list[x][[1]]) > 0) {
        # for each class comparison
        message(glue("Now plotting: {x}"))
        celltypeObj <- SubsetData(object,
                                  ident.use = x)
        Idents(celltypeObj) <- celltypeObj[[compare_by]]
        ipl <-
          lapply(
            names(combo_list[[x]]),
            function(y) {
              if (!is.null(y) & (!is.null(combo_list[[x]][[y]]))) {
                print(glue("Now plotting: {x} {y}"))
                bubbleplot(celltypeObj,
                           genes_plot = combo_list[[x]][[y]][, "gene"],
                           cluster_y = FALSE,
                           cluster_x = FALSE,
                           x_axis_title = str_wrap(y, width = 20),
                           y_axis_title = NULL,
                           do_return = TRUE,
                           pct_legend_title = "Percent\ngroup\nexpressing",
                           scale_legend_title = "Average\nscaled\nexpression",
                           ...)
              }
            }
          )
        grid.arrange(
          arrangeGrob(
            grobs = lapply(
              ipl,
              function(x) {
                x +
                  theme(
                    legend.position =
                      "none"
                  )
              }
            ),
            ncol = 1
          ),
          get_legend(ipl[[1]]),
          ncol = 2,
          widths = c(3, 1)
        )
      }
    }
  )
  return(opl)
}
