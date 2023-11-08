#' Generates pedigrees for grandmother example
#'
#' Optionally plots
#'
#' @param plot logical
#' @param pedNames character vector names of pedigrees
#' @param gm logical
#'
#' @return A list of pedigrees where `ped0` corresponds to H0: `unrelated`.
#'
#' @details Two examples are implemented, the grandmoter example (`gm = T`)
#' or with one pedigree added (`gm = F`).
#'
#' @examples
#' # peds = examplePedigrees(pedNames = NULL)
#' peds = examplePedigrees(pedNames = paste0("Ped", 0:2))
#' #peds = examplePedigrees(gm = F)
#' #peds = examplePedigrees(gm = F, pedNames = paste0("Ped", 0:3))
#' @export


examplePedigrees <- function(plot = T, pedNames = NULL, gm = T) {
  x <- swapSex(linearPed(2, sex = 1), 1, verbose = F)
  ped1 <- addSon(x, parents = c("3", "4"), id = "6")
  ped1 <- relabel(ped1, c("1", "EM1", "EM2", "EF1", "2", "3"))
  id <- c("1", "2", "3")
  x <- relabel(x, c(1, "EM1", "EM2", "EF2", "3"))
  ped2 <- addChildren(x, father = "EM2", id = "2", verbose = F)
  ped2 <- relabel(ped2, "EF1", "NN_1")
  ped2 <- relabel(ped2, c("2", "3"), c("3", "2"))
  x <- relabel(x, "EF1", "EF2")
  ped0 <- list(x, singleton("2", sex = 1))
  if(gm){
    peds <- list(ped0, ped1,  ped2)
    names(peds) = pedNames
    if(plot)
      plotPedList(peds, hatched = 1:3)
  }
  else{
    x = swapSex(linearPed(2, sex = 1), 1, verbose = F)
    x = relabel(x, c(1, "EM1", "EM2", "EF2", "3"))
    x = addParents(x, "1", "EM3", "EF3", verbose = F)
    y = swapSex(linearPed(3, sex = 1), 1, verbose = F)
    y = relabel(y,  c("EF3", paste(4:9)))
    ped3 = mergePed(x, y, by = "EF3")
    ped3 = relabel(ped3, "2", "9")
    ped3 = relabel(ped3, "EF1", "EF2")
    peds <- list(ped0, ped1, ped2, ped3)
    names(peds) = pedNames
    if(plot)
      plotPedList(peds, hatched = 1:3)
  }
  peds
}



