#' List the Containers Existed in Docker
#'
#' Utils function of getting the containers information in R.
#'
#' @param ... Nothing to input
#' @return A tibble
#' @export
#' @importFrom purrr map
#' @importFrom stats setNames
#' @importFrom utils read.delim
#'
#' @examples
#' # containers <- list_docker_container()
list_docker_container <- function(
    ...
){
  if(!requireNamespace("processx", quietly = TRUE)){
    message("Install processx...")
    install.packages('processx')
  }
  if(!requireNamespace("tibble", quietly = TRUE)){
    message("Install tibble...")
    install.packages('tibble')
  }
  columns <- c("ID", "Image", "Command", "Status", "Names")
  format <- paste(paste0("{{.", columns, "}}"), collapse = "\t")
  stdout <- processx::run("docker",
                          c("ps", "-all", paste0("--format=", format)),
                          echo = FALSE)$stdout

  if (stdout != "") {
    df <- utils::read.delim(text = stdout, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    df <- as.data.frame(df)
    colnames(df) <- columns
    df
  }else{
    purrr::map(columns, ~ character(0)) %>%
      stats::setNames(columns) %>%
      tibble::as_tibble()
  }
}
