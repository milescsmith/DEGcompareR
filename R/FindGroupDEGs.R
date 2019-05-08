#' @title FindGroupDEGs
#'
#' @description Find differentially expressed genes within an identity
#'   contrasted by a grouping variable (e.g. comparison between disease status
#'   or drug treatment (i.e. compare_by) within each cell type (i.e. ident_use).
#'   Performs 1-to-1 comparisons for each possible of combinations of the values
#'   of the grouping variable.
#'
#' @param object Processed Seurat scRNAseq object
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
#' @export
#' @importFrom gtools combinations
#' @importFrom glue glue
#' @importFrom Seurat SubsetData Idents Idents<- FindMarkers
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom purrr pluck map
#' @importFrom furrr future_map
#' @importFrom rlang enquo
FindGroupDEGs <- function(object,
                          ident_use,
                          compare_by,
                          test_use = "wilcox",
                          min_fold_change = 0,
                          min_pct_express_diff = 0,
                          pval_thresh = 0.05,
                          genes_of_interest = NULL,
                          cell_number_thresh = 10,
                          pval_use = p_val_adj) {
  p_val_adj <- NULL
  gene <- NULL
  avg_logFC <- NULL
  pct.1 <- NULL
  pct.2 <- NULL

  pval_use <- enquo(pval_use)
  if (is.null(genes_of_interest)) {
    genes_of_interest <- rownames(object)
  }
  combinations_results <-
    future_map(
      .x = unique(object[[ident_use]][[1]]),
      .progress = TRUE,
      .f = function(identity) {
        compareObj <- SubsetData(object = object,
                                 subset.name = ident_use,
                                 accept.value = identity)
        Idents(compareObj) <- compareObj[[compare_by]]

        unique_idents <- Idents(compareObj) %>% unique() %>% as.character()
        message(glue("Now processing {identity}"))

        if ((length(unique_idents) > 2)) {
          combos <- combinations(n = length(unique_idents),
                                 r = 2,
                                 v = unique_idents,
                                 repeats.allowed = FALSE)

          combo_DEs <- future_map(
            .x = 1:nrow(combos),
            .f = function(combo_row) {
              if (all(table(Idents(compareObj)) >= cell_number_thresh)) {
                message(glue("Now processing {identity} {combos[combo_row, 1]} {combos[combo_row, 2]}"))
                marks <- FindMarkers(object = compareObj,
                                     ident.1 = combos[combo_row, 1],
                                     ident.2 = combos[combo_row, 2],
                                     test.use = test_use) %>%
                  rownames_to_column("gene") %>%
                  filter((!!pval_use) < pval_thresh) %>%
                  filter(gene %in% genes_of_interest) %>%
                  filter(abs(avg_logFC) >= min_fold_change) %>%
                  filter(abs(pct.1 - pct.2) >= min_pct_express_diff)
                if (nrow(marks) > 0){
                  marks
                } else {
                  invisible()
                }

              }
            }
          )

          combo_names <- unlist(future_map(
            .x = 1:nrow(combos),
            .f = function(x) {
              if (all(table(Idents(compareObj)) >= cell_number_thresh)) {
                glue("{identity} {combos[x, 1]} vs {combos[x, 2]}")
              }
            }
          ))
          names(combo_DEs) <- combo_names
          combo_DEs <- combo_DEs[!sapply(combo_DEs, is.null)]
          combo_DEs
        } else if ((length(unique_idents) == 2) &
          all(table(Idents(compareObj)) >= cell_number_thresh)) {
          combo_DEs <- FindMarkers(object = compareObj,
                                   ident.1 = unique_idents[[1]],
                                   ident.2 = unique_idents[[2]],
                                   test.use = test_use) %>%
            rownames_to_column("gene") %>%
            filter((!!pval_use) < pval_thresh) %>%
            filter(gene %in% genes_of_interest) %>%
            filter(abs(avg_logFC) >= min_fold_change) %>%
            filter(abs(pct.1 - pct.2) >= min_pct_express_diff)
          combo_names <- glue("{identity} {unique_idents[[1]]} vs {unique_idents[[2]]}")
          names(combo_DEs) <- combo_names
          combo_DEs <- combo_DEs[!sapply(combo_DEs, is.null)]
          combo_DEs
        }
      }
    )

  message(object[[ident_use]] %>% unique() %>% pluck(1) %>% as.character())
  names(combinations_results) <- object[[ident_use]] %>% unique() %>% pluck(1) %>% as.character()
  combinations_results <- combinations_results[!sapply(combinations_results, is.null)]
  return(combinations_results)
}
