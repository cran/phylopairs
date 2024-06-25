#' Simulated Datasets
#'
#' @description
#' Simulated datasets are provided to help users understand functions. Lineage-pairs (n=190) were created from a simulated phylogenetic tree with 20 tips. A predictor variable was generated from random draws of a standard normal. Reponse variables of various structures were simulated as follows: 
#' 1. **Simulated Dataset 1** Unbounded response, linear relationship between response and predictor, no covariance in residuals
#' 2. **Simulated Dataset 2** Unbounded response, no relationship between response and predictor, no covariance in residuals
#' 3. **Simulated Dataset 3** Unbounded response, linear relationship between response and predictor, covariance in residuals
#' 4. **Simulated Dataset 4** Unbounded response, no relationship between response and predictor, covariance in residuals
#' 5. **Simulated Dataset 5** Bounded response, linear relationship between response and predictor on link scale, no covariance in residuals
#' 6. **Simulated Dataset 6** Bounded response, no relationship between response and predictor on link scale, no covariance in residuals
#' 7. **Simulated Dataset 7** Bounded response, linear relationship between response and predictor on link scale, covariance in residuals
#' 8. **Simulated Dataset 8** Bounded response, no relationship between response and predictor on link scale, covariance in residuals
#' 9. **Simulated Lineage-Pair Covariance Matrix** A 190*190 covariance matrix used in simulating the example datasets in this package.
#' 10. **Simulated Phylogenetic Tree** A 20-species ultrametric phylogenetic tree used in simulating the example datasets in this package.
#'
#' @format The datasets are as follows 
#' \describe{
#'   \item{data1}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data2}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data3}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data4}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data5}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data6}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data7}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{data8}{A data frame with 190 rows and 4 columns.
#'	   \describe{
#'       \item{sp1}{Identity of Species 1 in Pair}
#'       \item{sp2}{Identity of Species 2 in Pair}
#'       \item{pred}{Numeric predictor variable}
#'       \item{y}{Numeric response variable}
#'     }
#'   }
#'   \item{sim.cov.pairs}{A 190x190 covariance matrix
#'   }
#'   \item{sim.tree1}{A phylo object
#'   }
#'}
#'
#' @source Simulated data generated with the script provided in the 'inst' directory.
#'
#' @docType data
#' @name simulated.datasets
#' @aliases data1 data2 data3 data4 data5 data6 data7 data8 sim.cov.pairs sim.tree1
NULL

"data1"
"data2"
"data3"
"data4"
"data5"
"data6"
"data7"
"data8"
"sim.cov.pairs"
"sim.tree1"