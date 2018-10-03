#' @title FindGroupDEGs
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
FindGroupDEGs <- function(seuratObj,
                          ident_use,
                          compare_by,
                          test_use = "wilcox",
                          pval_thresh = 0.05,
                          genes_of_interest = NULL,
                          cell_number_thresh = 10,
                          pval_use = p_val_adj) {
  pval_use <- enquo(pval_use)
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

        if ((length(unique_idents) > 2)) {
          combos <- combinations(
            n = length(unique_idents),
            r = 2,
            v = unique_idents,
            repeats.allowed = FALSE
          )

          combo_DEs <- future_lapply(
            X = 1:nrow(combos),
            FUN = function(combo_row) {
              if (all(table(compareObj@ident) >= cell_number_thresh)) {
                print(glue("Now processing {identity} {combos[combo_row, 1]} {combos[combo_row, 2]}"))
                marks <- FindMarkers(
                  object = compareObj,
                  ident.1 = combos[combo_row, 1],
                  ident.2 = combos[combo_row, 2],
                  test.use = test_use
                ) %>%
                  rownames_to_column("gene") %>%
                  dplyr::filter((!!pval_use) < pval_thresh) %>%
                  dplyr::filter(gene %in% genes_of_interest)
                if (nrow(marks) > 0){
                  marks
                } else {
                  invisible()
                }

              }
            }
          )

          combo_names <- unlist(future_lapply(
            X = 1:nrow(combos),
            FUN = function(x) {
              if (all(table(compareObj@ident) >= cell_number_thresh)) {
                glue("{identity} {combos[x, 1]} vs {combos[x, 2]}")
              }
            }
          ))
          names(combo_DEs) <- combo_names
          combo_DEs <- combo_DEs[!sapply(combo_DEs, is.null)]
          combo_DEs
        } else if ((length(unique_idents) == 2) &
          all(table(compareObj@ident) >= cell_number_thresh)) {
          combo_DEs <- FindMarkers(
            object = compareObj,
            ident.1 = unique_idents[[1]],
            ident.2 = unique_idents[[2]],
            test.use = test.use
          ) %>%
            rownames_to_column("gene") %>%
            dplyr::filter((!!pval_use) < pval_thresh) %>%
            dplyr::filter(gene %in% genes_of_interest)
          combo_names <- glue("{identity} {unique_idents[[1]]} vs {unique_idents[[2]]}")
          names(combo_DEs) <- combo_names
          combo_DEs <- combo_DEs[!sapply(combo_DEs, is.null)]
          combo_DEs
        }
      }
    )
  names(combinations_results) <- unique(seuratObj@meta.data[, ident_use])
  combinations_results <- combinations_results[!sapply(combinations_results, is.null)]
  return(combinations_results)
}
