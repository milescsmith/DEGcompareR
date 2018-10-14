#' @title FindAllIntragroupMarkers
#'
#' @description Find differentially expressed genes within an identity
#'   contrasted by a grouping variable (e.g. comparison between disease status
#'   or drug treatment (i.e. compare_by) within each cell type (i.e. ident_use)
#'
#' @param seuratObj Processed Seurat scRNAseq object
#' @param ident_use Identity (from the ident slot or a meta.data column) by
#'   which to group the cells.
#' @param compare_by Meta.data column by which to group the cells within each
#'   identity for comparison.
#' @param min_fold_change Minimum log fold change difference between two groups
#' @param min_pct_express_diff Minimum difference in the percentage of each group
#'   that expresses the gene in question above a threshold
#' @param test_use Differential expression test to use.  Same values as in
#'   FindMarkers/FindAllMarkers (i.e wilcox, MAST, DESeq2, ...). Default:
#'   'wilcox'
#' @param pval_thresh Signifigance cutoff for DE testing. Default: 0.05
#' @param genes_of_interest A list of genes to use in DE testing. Default: NULL
#' @param cell_number_thresh Minimum number of cells that must be within an
#'   ident.use/group.by grouping before it is discarded. Default: 10
#' @param pval_use p value to use when filtering DE genes by the pval.thresh.
#'   Acceptable values are p_val_adj and p_val.  Default: p_val_adj
#'
#' @return Named list with genes
#'
#' @examples
#'
#' @export
#' @import magrittr
#' @importFrom gtools combinations
#' @importFrom glue glue
#' @importFrom future plan multiprocess
#' @importFrom future.apply future_lapply
#' @importFrom Seurat SubsetData SetAllIdent FindMarkers
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom rlang enquo
FindAllIntragroupMarkers <- function(seuratObj,
                                     ident_use,
                                     compare_by,
                                     test_use = "wilcox",
                                     min_fold_change = 0,
                                     min_pct_express_diff = 0,
                                     pval_thresh = 0.05,
                                     genes_of_interest = NULL,
                                     cell_number_thresh = 10,
                                     pval_use = p_val_adj) {
  plan(multiprocess)
  if (is.null(genes_of_interest)) {
    genes_of_interest <- rownames(seuratObj@raw.data)
  }
  combinations_results <-
    future_lapply(
      X = unique(seuratObj@meta.data[, ident_use]),
      FUN = function(identity) {
        compareObj <- SubsetData(
          object = seuratObj,
          subset.name = ident_use,
          accept.value = identity,
          subset.raw = TRUE
        )

        compareObj <- SetAllIdent(compareObj, compare_by)

        unique_idents <- as.character(unique(compareObj@ident))
        print(glue("Now processing {identity}"))



        combo_DEs <- FindAllMarkers(object = compareObj,
                                    genes.use = genes_of_interest,
                                    logfc.threshold = min_fold_change,
                                    test.use = test_use,
                                    min.diff.pct = min_pct_express_diff,
                                    min.cells.group = cell_number_thresh,,
                                    print.bar = TRUE,
                                    only.pos = TRUE,
                                    return.thresh = pval_thresh)
        combo_DEs <- combo_DEs[!sapply(combo_DEs, is.null)]
        combo_DEs
      })
  names(combinations_results) <- unique(seuratObj@meta.data[, ident_use])
  combinations_results <- combinations_results[!sapply(combinations_results, is.null)]
  return(combinations_results)
}
