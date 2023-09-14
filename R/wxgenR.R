#' \code{wxgenR} package
#'
#' A weather generator with seasonality
#'
#'
#' @docType package
#' @name wxgenR
#'
#' @rawNamespace import(stats, except = filter)
#' @import lubridate
#' @importFrom dplyr group_by summarise left_join glimpse mutate relocate if_else filter
#' @import msm
#' @import parallel
#' @import doParallel
#' @import foreach
#' @import sm
#' @import utils
#' @importFrom plyr ddply
#' @import magrittr
#'
#'
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c("temp", "state"))
