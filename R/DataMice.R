#' @name DataMice 
#' @title Mice Protein Expression Data Set
#'
#' @description The data set consists of the expression levels of 68 proteins that produced detectable signal in the nuclear fraction of 
#' cortex for a sample of 72 mice. There are 38 control mice and 34 trisomic mice. Several measurements were recorded for each protein
#' and for each mouse. The measurements containing missing observations in the original data were suppressed, so that one has between
#' 12 and 15 measurements per protein and per mouse. 
#' 
#' Mice may be further described based on the treatment they received (injected with memantine or saline), and on their behaviour 
#' (stimulated to learn or not). 
#'
#' @format A data frame of 72 rows (mice) and 905 columns (variables):
#' \itemize{
#'   \item{Protein_X_meas_Y}{numerical. The expression level for protein X at measurement Y. X has values between 1 and 68, Y has values
#'   between 1 and 12 or 15, according to the number of measurements.}
#'   \item{Genotype}{categorical. Two values: "Control" and "Ts65Dn" (trisomic mouse).}
#'   \item{Treatment}{categorical. Two values: "Memantine" and "Saline".}
#'   \item{Behaviour}{categorical. Two values: "C/S" (stimulated to learn) and "S/C" (not stimulated to learn).}
#'   \item{Class.mouse}{categorical. This variables creates eight classes of mice, based on crossing the categories of \code{Genotype},
#'   \code{Behaviour} and \code{Treatment}.}
#'   \item{MouseID}{factor. The key variable identifying each mouse in the sample.}
#' }
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression}
#' @references C. Higuera, K.J. Gardiner, and K.J. Cios (2015) Self-organizing feature maps identify proteins critical 
#' to learning in a mouse model of Down syndrome. PLoS ONE 10(6): e0129126.
"DataMice"

