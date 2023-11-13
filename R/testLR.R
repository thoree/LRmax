#' Finds classical likelihood ratio test statistics
#'
#' n named pedigrees are given. For one of this we have simulated marker data.
#' The test statistic Delta = -ln(R) is calculated. Here is the ratio of the
#' maximum likelihood (or LR) for pedigrees consistent
#' with H0 to to the overall maximum.
#' 1 to n,..., (n-1) to n and the maximum Z = max(LR.1.n, ..., LR.n-1.n)
#' are calculated
#'
#' @param pedigrees a list of pedigrees
#' @param sim pedigree with simulated marker data
#' @param truePed index or name of pedigree with marker data
#' @param ref index or name of pedigree in the denominator of each LR
#' @param H0 index or name of pedigree consistent with H0
#'
#' @return A data frame with LRs, the max and log likelihood of the first pedigree,
#' the reference and optionally a plot.
#'
#' @examples
#' pedigrees = examplePedigrees(plot = F, gm = F)
#' names(pedigrees) = paste0("Ped", 0:3)
#' truePed  = "Ped3"
#' ref = "Ped0"
#' nsim = 100
#' conditional = F
#' pedigrees[[truePed]] = setMarkers(pedigrees[[truePed]],
#'                   locusAttributes = NorwegianFrequencies[1:22])
#' pedigrees[[truePed]] = setMutmod(pedigrees[[truePed]], model = "proportional", rate = 0.001)
#' sim = simPedigree(pedigrees[[truePed]], idTarget ="2", idReferences = c("1", "3"),
#'                   nsim = nsim, seed = 1729, conditional = conditional)
#' H0 = c(1,4)
#' H0 = c("Ped0", "Ped3")
#' res = testLR(pedigrees, sim, truePed,  ref, H0)
#' res
#'
#' @export

testLR = function(pedigrees, sim, truePed = NULL, ref = 1, H0 = 1) {
  res = findLR(pedigrees, sim, truePed, ref)
  npeds = length(pedigrees)
  res2 = res[, 1:npeds]
  if(!is.integer(H0))
    H0 = (1:npeds)[names(pedigrees) %in% H0]
  if(length(H0) == 1)
    numerator = max(res2[,H0])
  else
    numerator = apply(res2[,H0], 1, max)
  denominator = apply(res2, 1, max)
  Delta = -2*log(numerator/denominator)
  Delta
}



