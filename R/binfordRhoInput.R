#' @name binfordRhoInput
#' @title An example input file for binfordRho
#' @description This dataset demonstrates the required input format for \code{binfordRho} when intervals are subsampled prior to analysis
#' The \code{binfordRhoInput} data frame contains the following columns: 
#' \itemize{
#'   \item INTTOP: depth (cm) to the top of the analyzed interval
#'   \item INTBOT: depth (cm) to the bottom of the analyzed interval
#'   \item WT_SED+BAG: mass (g) of the wet sediment and the container (i.e. Whirlpak bag)
#'   \item WT_VIAL: mass (g) of the empty vial used for freeze-drying
#'   \item WT_VIAL+SED: mass (g) of the vial after adding the subsample of wet sediment
#'   \item WT_FD: mass (g) of the vial plus sediment after freeze drying
#'   \item WT_GAMMA: mass (g) of freeze-dried sediment added to the gamma dating tube
#'   }
#'Note:
#'Columns 1-3 will have values for the entire core length
#'Columns 4-7 will only have values for the intervals prepared for dating by freeze-drying
#' @docType data
#' @usage data(binfordRhoInput)
#' @format a data frame containing 60 observations of 7 variables
#' @author Adam Jeziorski
NULL
"binfordRhoInput"
