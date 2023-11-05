#' Simulation for grandmother case
#'
#' Conditional and unconditional simulation for hard coded example
#'
#' @param pedigrees a list of three pedigrees
#' @param Nsim integer
#' @param seed integer
#' @param conditional logical
#' @param character `Ped1`, `Ped2`, or `Ped3`
#'
#' @return A data frame with LR.1.3, LR.2.3 and max, i.e., LRmax
#'
#' @examples
#' Set1 <- c(
#'   "CSF1PO", "D2S1338", "D3S1358", "D5S818", "D7S820",
#'   "D8S1179", "D13S317", "D16S539", "D18S51", "D19S433", "D21S11", "FGA",
#'   "TH01", "TPOX", "VWA"
#' )
#' db <- forrel::NorwegianFrequencies[Set1]
#' nsim = 2
#' seed17 = 17
#' cond = F
#' mut = T#'
#'
#' # Generate pedigrees
#' peds = gmPedigrees(plot = F)
#'
#' # Simulate from Ped1
#' peds[[1]] <- setMarkers(peds[[1]], locusAttributes = db)
#' if(mut) peds[[1]] <- setMutmod(peds[[1]], model = "proportional", rate = 0.001)
#' LR1 <- gmExample(
#'   pedigrees = peds, Nsim = nsim, seed = seed17,
#'   conditional = cond, simulateFrom = "Ped1"
#' )
#' log1Max <- log10(LR1[, 3])
#'
#' # Simulate from Ped2
#' peds[[2]] <- setMarkers(peds[[2]], locusAttributes = db)
#' if(mut) peds[[2]] <- setMutmod(peds[[2]], model = "proportional", rate = 0.001)
#' LR2 <- gmExample(
#'   pedigrees = peds, Nsim = nsim, seed = seed17,
#'   conditional = cond, simulateFrom = "Ped2"
#' )
#' log2Max <- log10(LR2[, 3])
#'
#' # Simulate from Ped3
#' peds[[3]] <- setMarkers(peds[[3]], locusAttributes = db)
#' if (mut) peds[[3]] <- setMutmod(peds[[3]], model = "proportional", rate = 0.001)
#' LR3 <- gmExample(
#'   pedigrees = peds, Nsim = nsim, seed = seed17,
#'   conditional = cond, simulateFrom = "Ped3"
#' )
#' log3Max <- log10(LR3[, 3])

#' boxplot(log1Max, log2Max, log3Max, names = c("Ped1", "Ped2", "Ped3"),  ylab = "log10(LRmax)")
#' criticalValue = quantile(log3Max, probs = 0.95)
#' abline(h = criticalValue)
#' @export


gmExample <- function(pedigrees, Nsim, seed, conditional = F, simulateFrom = "Ped3") {
  set.seed(seed)
  if (conditional) {
    sim <- profileSim(pedigrees[[simulateFrom]], N = 1, ids = c("1", "3")) |>
      profileSim(N = Nsim, ids = "2")
  } else {
    sim <- profileSim(pedigrees[[simulateFrom]], N = Nsim, ids = c("1", "2", "3"))
  }


  LRs <- matrix(nrow = nsim, ncol = 3)
  if (simulateFrom == "Ped1") {
    for (i in 1:Nsim) {
      LRs[i, 1:2] <- as.double(kinshipLR(
        "Ped1" = sim[[i]],
        "Ped2" = pedigrees[[2]],
        "Ped3" = pedigrees[[3]], ref = 3,
        source = simulateFrom
      )$LRtotal[1:2])
    }
  } else if (simulateFrom == "Ped2") {
    for (i in 1:Nsim) {
      LRs[i, 1:2] <- as.double(kinshipLR(
        "Ped1" = pedigrees[[1]],
        "Ped2" = sim[[i]],
        "Ped3" = pedigrees[[3]], ref = 3,
        source = simulateFrom
      )$LRtotal[1:2])
    }
  } else if (simulateFrom == "Ped3") {
    for (i in 1:Nsim) {
      LRs[i, 1:2] <- as.double(kinshipLR(
        "Ped1" = pedigrees[[1]],
        "Ped2" = pedigrees[[2]],
        "Ped3" = sim[[i]], ref = 3,
        source = simulateFrom
      )$LRtotal[1:2])
    }
  }
  LRmax <- pmax(LRs[, 1], LRs[, 2])
  res <- data.frame(LR.1.3 = LRs[, 1], LR.2.3 = LRs[, 2], LRmax = LRmax)
  res
}



