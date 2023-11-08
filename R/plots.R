#' Plot LRs and max
#'
#' Based on output from `findLR()` a boxplot of the LRs and Z or a density plot of
#' is produced
#'
#' @param res a data frame from`findLR()`
#' @param plot character, `box`, `density` or NULL
#' @param ref integer or character indicating reference pedigree
#' @param ref integer or character indicating true ped
#'
#' @return Various plots
#'
#' @examples
#' peds = examplePedigrees(plot = F, gm = F)
#' names(peds) = paste0("Ped", 0:3)
#' truePed  = "Ped0"
#' ref = "Ped0"
#' peds[[truePed]] = setMarkers(peds[[truePed]],
#'                   locusAttributes = NorwegianFrequencies[1:22])
#' peds[[truePed]] = setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = simPedigree(peds[[truePed]], idTarget ="2", idReferences = c("1", "3"),
#'                   nsim = 10, seed = 1729, conditional = F)
#'
#' res = findLR(peds, sim,  truePed = truePed,  ref = ref)
#' plots(res, plot = "box", ref = ref)
#' plots(res, plot = "density", ref = ref)
#' @export


plots =  function(res, plot = "box", ref = 1, truePed = 1, tit = NULL) {
  d = dim(res)
  npeds = d[2] - 2 # number of pedigrees
  cn = colnames(res)[1:(npeds+1)] #last is Z
  if(is.integer(ref))
    index = (1:npeds)[-ref]
  else
    index = (1:npeds)[names(peds) != ref]

  if(plot == "box"){
    boxplot(log10(res[,c(index, npeds + 1)]), names = c(cn[index], "Z"), ylab = "log10(LR.i.0)")
    title(tit)
  }
  if(plot == "density"){
    plot(density(log10(res$Z)), xlab = "log10(Z)", main = tit)
         # main = paste("Distr: log10(Z), ", nsim, " sims from ", truePed))
  }

}



