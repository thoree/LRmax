#' Performs max LR test
#'
#' n named pedigrees are given, one with marker data. The LRs comparing
#' 1 to n,..., (n-1) to n and the maximum Z = max(LR.1.n, ..., LR.n-1.n)
#' are calculated
#'
#' @param pedigrees a list of pedigrees
#' @param truePed integer, index of pedigree with marker data
#'
#' @return A data frame with LRs and the max
#'
#' @examples
#' peds = gmPedigrees(plot = F)
#' truePed  = 2
#' nsim =2
#' conditional = T
#' peds[[truePed]] = setMarkers(peds[[truePed]],
#'                   locusAttributes = NorwegianFrequencies[1:5])
#' peds[[truePed]] = setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = simPedigree(peds[[truePed]], idTarget ="2", idReferences = c("1", "3"),
#'                   nsim = nsim, seed = 1729, conditional = conditional)
#'
#' LRmaxTest(peds, sim, nsim = nsim, truePed = truePed)
#'
#' @export


LRmaxTest <- function(pedigrees, sim, nsim = NULL, truePed = NULL) {
  npeds = length(pedigrees)
  if(is.null(truePed))
    truePed = npeds
  res = matrix(nrow = nsim, ncol = npeds + 1)
  if(nsim == 1){
    pedigrees[[truePed]] = transferMarkers(sim, pedigrees[[truePed]])
    lrs = kinshipLR(pedigrees, ref = npeds, source = truePed)$LRtotal
    res[1, 1:npeds] = as.double(lrs)
  } else {
    for (i in 1:nsim){
      pedigrees[[truePed]] = transferMarkers(sim[[i]], pedigrees[[truePed]])
      lrs = kinshipLR(pedigrees, ref = npeds, source = truePed)$LRtotal
      res[i, 1:npeds] = as.double(lrs)
    }
  }
  if(nsim == 1)
    Z = max(res[1,1:(npeds-1)])
  else
    Z = apply(res[,1:(npeds-1)], 2, max)

  res[, npeds+1] = Z
  colnames(res) = c(names(lrs), "Z")
  res = data.frame(res)
  res
}



