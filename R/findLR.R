#' Finds LRs and the maximum
#'
#' n named pedigrees are given and simulated marker data. The LRs comparing
#' 1 to n,..., (n-1) to n and the maximum Z = max(LR.1.n, ..., LR.n-1.n)
#' are calculated
#'
#' @param pedigrees a list of pedigrees
#' @param sim pedigree with simulated marker data
#' @param truePed index or name of pedigree with marker data
#' @param ref index or name of pedigree in the numerator of LR
#'
#' @return A data frame with LRs, the max and log likelihood of the first pedigree,
#' the reference and optionally a plot.
#'
#' @examples
#' peds = examplePedigrees(plot = F, gm = F)
#' names(peds) = paste0("Ped", 0:3)
#' truePed  = "Ped0"
#' ref = "Ped1"
#' nsim = 2
#' conditional = F
#' peds[[truePed]] = setMarkers(peds[[truePed]],
#'                   locusAttributes = NorwegianFrequencies[1:22])
#' peds[[truePed]] = setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = simPedigree(peds[[truePed]], idTarget ="2", idReferences = c("1", "3"),
#'                   nsim = nsim, seed = 1729, conditional = conditional)
#'
#' res = findLR(peds, sim, truePed = truePed,  ref = ref)
#'
#' @export

findLR = function(pedigrees, sim, truePed = NULL, ref = 1) {
  npeds = length(pedigrees)
  nsim = length(sim)
  # Convert ref to integer if necessary
  if(!is.integer(ref))
    ref = (1:npeds)[names(pedigrees) == ref]

  if(is.null(truePed))
    truePed = npeds
  res = matrix(nrow = nsim, ncol = npeds + 2)
  if(nsim == 1){
    pedigrees[[truePed]] = transferMarkers(sim, pedigrees[[truePed]])
    lrs = kinshipLR(pedigrees, ref = ref, source = truePed)
    res[1, 1:npeds] = as.double(lrs$LRtotal)
    res[1, npeds+2] =  sum(log(lrs$likelihoodsPerMarker[,ref]))
  } else {
    for (i in 1:nsim){
      pedigrees[[truePed]] = transferMarkers(sim[[i]], pedigrees[[truePed]])
      lrs = kinshipLR(pedigrees, ref = ref, source = truePed)
      res[i, 1:npeds] = as.double(lrs$LRtotal)
      res[i, npeds+2] =  sum(log(lrs$likelihoodsPerMarker[,ref]))
    }
  }

  index = (1:npeds)[-ref]
  if(nsim == 1)
    Z = max(res[1,index])
  else
    Z = apply(res[, index], 1, max)

  res[, npeds+1] = Z
  cn = names(lrs$LRtotal)
  colnames(res) = c(cn, "Z", "logRef")
  res = data.frame(res)
  res
}



