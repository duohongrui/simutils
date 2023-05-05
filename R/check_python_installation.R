#' Check Python Installation
#'
#' @param ... NULL
#' @importFrom stringr str_split
#' @importFrom utils compareVersion
#' @export
#'
#' @examples
#' check_python_installation()
check_python_installation <- function(...){
  if(!requireNamespace("crayon", quietly = TRUE)){
    stop("Package \"crayon\" must be installed by \"install.packages('crayon')\" command.")
  }
  if(!requireNamespace("processx", quietly = TRUE)){
    utils::install.packages("processx")
  }
  ### check python version
  python_symbol <- "python3"
  version_out <- processx::run(python_symbol,
                               c("--version"),
                               error_on_status = FALSE,
                               timeout = 2)
  if(version_out$status != 0 & version_out$stdout == ""){
    python_symbol <- "python"
    version_out <- processx::run(python_symbol,
                                 c("--version"),
                                 error_on_status = FALSE,
                                 timeout = 2)
    if(version_out$status != 0 & version_out$stdout == ""){
      stop(crayon::red(paste0("\u274C Python is not installed on your device.")))
    }
  }
  message(crayon::green(paste0("\u2714 Python is already installed.")))
  # version
  version <- stringr::str_split(version_out$stdout, pattern = " ", simplify = T)[2] %>% trimws()
  if(utils::compareVersion("3.6.0", version) > 0){
    stop(crayon::red(paste0("\u274C Your python is ", version, ". But 3.6.0 or higher is required.")))
  }else{
    message(crayon::green(paste0("\u2714 Your python version is satisfied.")))
  }

  ### check python modules
  if(python_symbol == "python"){
    pip_com <- "pip"
  }else{
    pip_com <- "pip3"
  }

  depend_packages <- c("numpy",
                       "scipy",
                       "pandas",
                       "newick",
                       "prosstt")
  for(pa in depend_packages){
    modules_out <- processx::run(pip_com,
                                 c("show", pa),
                                 error_on_status = FALSE,
                                 timeout = 10)
    if(modules_out$status == 0){
      message(crayon::green(paste0("\u2714 ", pa, " module is installed.")))
    }else{
      if(pa == "prosstt"){
        message(crayon::red(paste0("\u274C prosstt is not installed, please refer to https://github.com/soedinglab/prosstt for instructions.")))
      }else{
        message(crayon::red(paste0("\u274C ", pa, " module is not installed")))
      }
    }
  }

}

print_error <- function(x, proc){
  print(x)
}
