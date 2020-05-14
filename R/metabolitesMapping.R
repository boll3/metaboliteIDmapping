###########################################################
#' ID Mapping table of nine different metabolite ID formats
#'
#' Four different sources of annotated metabolites, i.e., \code{HMDB},
#' \code{ChEBI}, \code{CompTox}, and the \code{graphite} R package, 
#' have been retrieved to compile a comprehensive mapping of available
#' metabolite IDs. ID formats that are represented in the mapping table are:
#' DTXCID (Comptox),
#' DTXSID (Comptox),
#' CAS-number,
#' CID (Pubchem),
#' SID (Pubchem),
#' HMDB,
#' ChEBI,
#' KEGG,
#' and Drugbank
#'
#' @name metabolitesMapping
#'
#' @docType data
#'
#' @usage metabolitesMapping
#'
#' @format A tibble with 9 variables and over 1.1 million metabolites:
#' \describe{
#'     \item{DTXCID}{DSSTox structure identifier, character}
#'     \item{DTXSID}{DSSTox substance identifier, character}
#'     \item{CAS}{CAS registry number, character}
#'     \item{CID}{Pubchem compound identifier, character}
#'     \item{CID}{Pubchem substance identifier, character}
#'     \item{HMDB}{Human Metabolome Database identifier (new format), character}
#'     \item{ChEBI}{Chemical Entities of Biological Interest identifier, character}
#'     \item{KEGG}{KEGG Compound identifier, character}
#'     \item{Drugbank}{Drugbank identifier, character}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' metabolitesMapping
NULL
###########################################################