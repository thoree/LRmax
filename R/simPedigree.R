#' Simulates marker data for a pedigree
#'
#' Unconditional or conditional simulation.
#'
#' @param ped a pedigree or list of pedigrees with markers attached
#' @param idTarget character POI
#' @param idReferences character vector References
#' @param integer nsim
#' @param integer seed
#' @param logic conditional

#'
#' @return Pedigree with simulations
#'
#' @examples
#' peds = gmPedigrees(plot = F)
#' truePed = 3
#' peds[[truePed]] <- setMarkers(peds[[truePed]],
#'              locusAttributes =  NorwegianFrequencies)
#' peds[[truePed]] <- setMutmod(peds[[truePed]], model = "proportional", rate = 0.001)
#' sim = simPedigree(peds[[truePed]], idTarget = "2", idReferences = c("2", "3"),
#'               nsim = 2, conditional = T)
#' sim
#'
#' @export


simPedigree <- function(ped, idTarget, idReferences, nsim = 1, seed = NULL,
                        conditional = F) {
  if (conditional) {
    sim <- profileSim(ped, N = 1, ids = idReferences) |>
      profileSim(N = nsim, ids = idTarget, seed = seed)
  } else {
    sim <- profileSim(ped, N = nsim,
                      ids = c(idTarget, idReferences), seed = seed)
  }

}



