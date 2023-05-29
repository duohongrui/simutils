
<img src="man/figures/simutils_logo.png" align="right" width = "167px" height="193px"/>

# Simutils implements essential utils functions for estimating parameters from real data, simulating and evaluating the simulated single-cell RNA-seq datastets.

Simutils contains many useful utils functions to establish the standard
pipeline for defining a simulation method, preprocessing reference
datasets with specific requirements, preparing important files for
special simulators, evaluating the performance of different simulation
methods and so on.

## Installation

You can install the development version of simutils from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("duohongrui/simutils")
```

``` r
library(simutils)
```

    ## Registered S3 method overwritten by 'SeuratDisk':
    ##   method            from  
    ##   as.sparse.H5Group Seurat

    ## Welcome to simutils

## Check Python Installation

Simutils contains a function for checking the Python installation and
environmrnt configuration for
[**PROSSTT**](https://github.com/duohongrui/simmethods/blob/master/R/28-PROSSTT.R)
method.

``` r
simutils::check_python_installation()
```

    ## ✔ Python is already installed.

    ## ✔ Your python version is satisfied.

    ## ✔ numpy module is installed.

    ## ✔ scipy module is installed.

    ## ✔ pandas module is installed.

    ## ✔ newick module is installed.

    ## ✔ prosstt module is installed.

## Match Cells From Real And Simulated Data

We adopted the Hungarian algorithm to match the cells from reference and
simulated datasets. In addition, we also provide an improved Hungarian
algorithm.

``` r
set.seed(666)
ref_data <- matrix(rpois(10 ^ 6, 2),
                   ncol = 1000,
                   nrow = 1000,
                   dimnames = list(paste0("ref_gene", 1:1000),
                                   paste0("ref_cell", 1:1000)))
set.seed(666)
sim_data <- matrix(rpois(10 ^ 6, 2.5),
                   ncol = 1000,
                   nrow = 1000,
                   dimnames = list(paste0("sim_gene", 1:1000),
                                   paste0("sim_cell", 1:1000)))
```

``` r
match_result <- simutils::match_cells(ref_data = ref_data,
                                      sim_data = sim_data,
                                      t = TRUE,
                                      algorithm = "Hungarian")
```

    ## Performing PCA...
    ## Performing Harmony...
    ## Calculate correlation matrix...

    ## 

    ## Match simulated and real cells using Hungarian...
    ##     reference simulation match_value
    ## 1  ref_cell77  sim_cell1   0.3809364
    ## 2 ref_cell777  sim_cell2   0.4301080
    ## 3  ref_cell45  sim_cell3   0.5085714
    ## 4 ref_cell732  sim_cell4   0.4225210
    ## 5 ref_cell353  sim_cell5   0.4255942
    ## 6 ref_cell818  sim_cell6   0.4712125

``` r
head(match_result[["cell_pair"]][order(match_result[["cell_pair"]]$match_value, decreasing = TRUE), ])
```

    ##       reference  simulation match_value
    ## 921 ref_cell487 sim_cell921   0.6337095
    ## 20  ref_cell213  sim_cell20   0.6115246
    ## 186 ref_cell775 sim_cell186   0.5933733
    ## 650 ref_cell495 sim_cell650   0.5933733
    ## 70   ref_cell62  sim_cell70   0.5922209
    ## 840 ref_cell660 sim_cell840   0.5893397

We can also use improved Hungarian algorithm:

``` r
match_result2 <- simutils::match_cells(ref_data = ref_data,
                                       sim_data = sim_data,
                                       t = TRUE,
                                       algorithm = "Improved_Hungarian")
```

    ## Performing PCA...
    ## Performing Harmony...
    ## Calculate correlation matrix...
    ## Match simulated and real cells using improved Hungarian...
    ##   reference  simulation match_value
    ## 1 ref_cell1 sim_cell427   0.4751501
    ## 2 ref_cell2 sim_cell851   0.5117407
    ## 3 ref_cell3 sim_cell125   0.5028091
    ## 4 ref_cell4 sim_cell935   0.4354862
    ## 5 ref_cell5 sim_cell111   0.4327971
    ## 6 ref_cell6 sim_cell671   0.4223289

``` r
head(match_result2[["cell_pair"]][order(match_result2[["cell_pair"]]$match_value, decreasing = TRUE), ])
```

    ##       reference  simulation match_value
    ## 487 ref_cell487 sim_cell921   0.6337095
    ## 213 ref_cell213  sim_cell20   0.6115246
    ## 495 ref_cell495 sim_cell650   0.5933733
    ## 775 ref_cell775 sim_cell186   0.5933733
    ## 62   ref_cell62  sim_cell70   0.5922209
    ## 660 ref_cell660 sim_cell840   0.5893397
