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
#' Set1 <- c(
#'   "CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820",
#'   "D8S1179", "D13S317", "D16S539", "D18S51", "D19S433", "D21S11", "FGA",
#'   "TH01", "TPOX", "VWA"
#' )
#' db <- forrel::NorwegianFrequencies[Set1]
#' mut = T
#'
#' # Simulate once from Ped truePed
#' nsim = 1
#' truePed = 1
#' peds[[truePed]] <- setMarkers(peds[[truePed]], locusAttributes = db)
#' if (mut) peds[[truePed]] <- setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = profileSim(peds[[truePed]], N = 1, ids = c("1", "2", "3"), seed = 1729)
#' pedigrees = peds
#' LRmaxTest(peds, sim, nsim = nsim, truePed = truePed)
#'
#' nsim = 2
#'
#' peds[[truePed]] <- setMarkers(peds[[truePed]], locusAttributes = db)
#' if (mut) peds[[truePed]] <- setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = profileSim(peds[[truePed]], N = nsim, ids = c("1", "2", "3"), seed = 1729)
#' #pedigrees = peds
#' LRmaxTest(peds, sim, nsim = nsim)
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
    cn = c(names(lrs), "Z")
    lrs = as.double(lrs)
    Z = max(lrs[-npeds])
    res[1, ] = c(lrs, Z)
    dimnames(res) = list("obs", cn)
    res = data.frame(res)
  } else {
    for (i in 1:nsim){
      pedigrees[[truePed]] = transferMarkers(sim[[i]], pedigrees[[truePed]])
      lrs = kinshipLR(pedigrees, ref = npeds, source = truePed)$LRtotal
      lrs2 = as.double(lrs)
      Z = max(lrs2[-npeds])
      res[i, ] = c(lrs2, Z)
    }
    cn = c(names(lrs), "Z")
    dimnames(res) = list(1:nsim, cn)
    res = data.frame(res)
   }
  res
}



