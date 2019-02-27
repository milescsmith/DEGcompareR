#' WriteGroupDEGs
#'
#' Given a list from FindGroupDEGs, write the expression values and percentage
#' of cells expressing to csv files
#'
#' @param combo_list Differentially expressed gene list for combination of contrast factors from FindGroupDEGs
#' @param location Location to which to write the csv files
#'
#' @return
#' @export
#'
#' @importFrom future.apply future_lapply
#' @importFrom stringr str_sub str_replace_all
#' @importFrom glue glue
#' @importFrom data.table fwrite
#'
#' @examples
WriteFindGroupDEGs <- function(combo_list, location) {
  if (!str_sub(location, -1) == "/") {
    location <- glue("{location}/")
  }
  # for each cell type
  future_lapply(
    names(combo_list),
    function(x) {
      future_lapply(
        names(combo_list[[x]]),
        function(y) {
          if (!is.null(y) & (!is.null(combo_list[[x]][[y]]))) {
            y_name <- str_replace_all(y,
              pattern = " ",
              replacement = "_"
            )
            y_name <- str_replace_all(y_name,
              pattern = "\\/",
              replacement = "_"
            )
            y_name <- str_replace_all(y_name,
              pattern = "\\+",
              replacement = "-pos"
            )
            y_name <- str_replace_all(y_name,
              pattern = "\\-",
              replacement = "-neg"
            )
            fwrite(
              x = combo_list[[x]][[y]],
              file = glue("{location}{y_name}")
            )
          }
        }
      )
    }
  )
}
