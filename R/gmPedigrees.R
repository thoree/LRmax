#' Generates pedigrees for grandmother example
#'
#' Optionally plots
#'
#' @param plot logical
#' @param tit title Character of length 3
#'
#' @return A list of pedigrees
#'
#' @examples
#' peds = gmPedigrees()
#'
#' @export


gmPedigrees <- function(plot = T, tit = c("Ped1", "Ped2", "Ped3")) {
  x <- swapSex(linearPed(2, sex = 1), 1, verbose = F)
  ped1 <- addSon(x, parents = c("3", "4"), id = "6")
  ped1 <- relabel(ped1, c("1", "EM1", "EM2", "EF1", "2", "3"))
  id <- c("1", "2", "3")
  x <- relabel(x, c(1, "EM1", "EM2", "EF2", "3"))
  ped2 <- addChildren(x, father = "EM2", id = "2", verbose = F)
  ped2 <- relabel(ped2, "EF1", "NN_1")
  ped2 <- relabel(ped2, c("2", "3"), c("3", "2"))
  x <- relabel(x, "EF1", "EF2")
  ped3 <- list(x, singleton("2", sex = 1))
  peds <- list("Ped1" = ped1, "Ped2" = ped2, "Ped3" = ped3)
  if(plot)
    plotPedList(peds, hatched = 1:3, titles = tit)
  peds
}



