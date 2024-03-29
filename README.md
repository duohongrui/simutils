
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

Here we demonstrate some useful functions in simutils package:

<ul>
<li>
<a href="#a1">Check Python Installation</a>
</li>
<li>
<a href="#a2">Match Cells From Real And Simulated Data</a>
</li>
<li>
<a href="#a3">Format Conversion of Single-Cell Data</a>
</li>
<li>
<a href="#a4">Make Lineage Trees for Single-Cell Data</a>
</li>
</ul>
<h2>
<a name="a1">Check Python Installation</a>
</h2>

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

<h2>
<a name="a2">Match Cells From Real And Simulated Data</a>
</h2>

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

<h2>
<a name="a3">Format Conversion of Single-Cell Data</a>
</h2>

Simsite provides the function of converting single-cell data formats
from *SingleCellExperimental* to *Seurat*, *list* and *h5ad*.

Here we construct a data with *SingleCellExperimental* format:

``` r
library(SingleCellExperiment)
```

``` r
set.seed(111)
data <- matrix(rpois(10 ^ 6, 2),
               ncol = 1000,
               nrow = 1000,
               dimnames = list(paste0("ref_gene", 1:1000),
                               paste0("ref_cell", 1:1000)))
SCE <- SingleCellExperiment::SingleCellExperiment(list(counts = data),
                                                  colData = data.frame("cell_name" = colnames(data)),
                                                  rowData = data.frame("gene_name" = rownames(data)))
SCE
```

    ## class: SingleCellExperiment 
    ## dim: 1000 1000 
    ## metadata(0):
    ## assays(1): counts
    ## rownames(1000): ref_gene1 ref_gene2 ... ref_gene999 ref_gene1000
    ## rowData names(1): gene_name
    ## colnames(1000): ref_cell1 ref_cell2 ... ref_cell999 ref_cell1000
    ## colData names(1): cell_name
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

To Seurat object:

``` r
Seurat <- simutils::data_conversion(SCE_object = SCE, return_format = "Seurat")
Seurat
```

    ## An object of class Seurat 
    ## 1000 features across 1000 samples within 1 assay 
    ## Active assay: originalexp (1000 features, 0 variable features)

To h5ad file (and you will get a path at which the file locates):

``` r
h5ad <- simutils::data_conversion(SCE_object = SCE, return_format = "h5ad")
```

    ## Creating h5Seurat file for version 3.1.5.9900

    ## Adding counts for originalexp

    ## Adding data for originalexp

    ## No variable features found for originalexp

    ## Adding feature-level metadata for originalexp

    ## Validating h5Seurat file

    ## Adding data from originalexp as X

    ## Transfering meta.features to var

    ## Adding counts from originalexp as raw

    ## Transfering meta.features to raw/var

    ## Transfering meta.data to obs

    ## Your data has been save to /var/folders/1l/xmc98tgx0m37wxtbtwnl6h7c0000gn/T//Rtmp9kqkJv/20230529193826.h5ad

``` r
h5ad
```

    ## $file_type
    ## [1] "h5ad"
    ## 
    ## $save_path
    ## [1] "/var/folders/1l/xmc98tgx0m37wxtbtwnl6h7c0000gn/T//Rtmp9kqkJv/20230529193826.h5ad"

<h2>
<a name="a4">Make Lineage Trees for Single-Cell Data</a>
</h2>

The tree lineage of cell clusters within single-cell data can be traced
by `make_trees` function in simutils.

``` r
set.seed(555)
data <- cbind(matrix(rpois(2.5e5, 2), ncol = 250, nrow = 1000),
              matrix(rpois(2.5e5, 4), ncol = 250, nrow = 1000),
              matrix(rpois(2.5e5, 8), ncol = 250, nrow = 1000),
              matrix(rpois(2.5e5, 12), ncol = 250, nrow = 1000))
colnames(data) <- paste0("Cell", 1:ncol(data))
rownames(data) <- paste0("Gene", 1:nrow(data))
```

Newick format:

``` r
newick <- simutils::make_trees(ref_data = data,
                               is_Newick = TRUE,
                               is_parenthetic = FALSE,
                               return_group = TRUE)
```

    ## Loading required package: amap

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Your data has 3 groups

``` r
newick$phyla
```

    ## [1] "(group2:3.5470224281145,(group1:5.70650677442804,group3:5.70650677442804):3.5470224281145);"

Phylo format:

``` r
phylo <- simutils::make_trees(ref_data = data,
                              is_Newick = FALSE,
                              is_parenthetic = TRUE,
                              return_group = TRUE)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

    ## Your data has 3 groups

``` r
str(phylo$phyla)
```

    ## List of 1
    ##  $ :List of 4
    ##   ..$ edge       : int [1:4, 1:2] 4 4 5 5 1 5 2 3
    ##   ..$ edge.length: num [1:4] 3.55 3.55 5.71 5.71
    ##   ..$ Nnode      : int 2
    ##   ..$ tip.label  : chr [1:3] "group2" "group1" "group3"
    ##   ..- attr(*, "class")= chr "phylo"
    ##   ..- attr(*, "order")= chr "cladewise"
