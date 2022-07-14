#' Local beta coefficients from source sites
#'
#' A 25*15 numeric matrix containing the local estimated beta coefficients from 15 source sites.
#' Rows stand for different variables while columns mean different source sites.
#'
"B_source"

#' Local beta coefficients from source sites
#'
#' A 25*16 numeric matrix containing the local estimated beta coefficients from all 16 source sites.
#' Each would become a target source under the federated learning framework. Rows stand for different
#' variables while columns mean different sites.
#'
"B_all"


#' Data covariance matrix from the target site.
#'
#' A 25*25 covariance matrix calculated from the target site.
#'
"Sigma_target"

#' Data covariance matrix from each site.
#'
#' A list containing 16 elements, each of which being a 25*25 covariance
#' matrix calculated from one site.
#'
"Sigma_all"


#' Time to event for the validation dataset.
#'
#' A numeric vector denoting the followup time for 2000 patients in the validation dataset.
#'
"x.valid"

#' Baseline covariates in the validation dataset.
#'
#' A numeric matrix including 25 baseline covariates for the validation dataset.
#'
"z.valid"

#' Status indicator for the validation dataset.
#'
#' A 0-1 vector indicating the status of 2000 patients in the validation dataset. 1=dead, 0=censored.
#'
"delta.valid"
