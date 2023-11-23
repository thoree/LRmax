#' Finds Generalised Likelihood Ratios (GLRs)
#'
#' Based on a `dviData` objects, the GLR for specified tests are computed
#'
#' @param dvi A `dviData` object.
#' @param pairings List. See details.
#' @param list Hypotheses.
#' @param disableMutations A logical, or NA (default). The default action is to
#'   disable mutations in all reference families without Mendelian errors.
#' @param verbose A logical.
#'
#' @return A data frame with hypotheses tested, max loglikelihid in numerator and
#'   denominator and GLRs
#'
#' @details The Generalised Likelihood Ratio (GLR) statistic is defined as the
#'   ratio of the maximum likelihood for the alternatives in the numerator
#'   to the maximum in the denominator. The default`dvir::generatePairings(dvi)`
#'   gives all hypotheses. Specific tests can be specified as shown in an example:
#'   `pairings = list(V1 = "M1")` gives a test for H0: V1 = M1 against H1: V1 != M1.
#'
#' @examples
#' dviGLR(example2, pairings = list(V1 = "M1"))
#'
#' # All tests with output from dviJoint
#' dviGLR(example2, verbose = T)
#'
#' @export

dviGLR = function(dvi, pairings = generatePairings(dvi),
                  disableMutations = FALSE, verbose = FALSE){
  if(class(dvi) != "dviData")
    stop("First argument needs to be a dviData object")

  # Do joint analysis and find max loglikelihood
  r = dviJoint(dvi, verbose = F, disableMutations = disableMutations)
  if(verbose)
    print(r)
  loglik = r$loglik
  na = dim(r)[1]

  # The number of hypotheses to be tested
  nhyp = sum(unlist(lapply(pairings, function(x) length(x))))

  # Initialise log likelihoods corresponding to numerator and denominator
  l0 = l1 = rep(-Inf, nhyp)

  # To contain names of null hypotheses in output:
  hyp = rep(NA, nhyp)

  # Loop through hypotheses for each victim
  npairs = length(pairings)
  s = 0
  for (i in 1:npairs){
    victim = names(pairings[i])
    for(j in pairings[[i]]){
      s = s +1
      hyp[s] = paste(victim, j, sep = " = ")

      # Find indices for numerator  loglikelihoods
      I0 = (1:na)[r[,victim] == j]
      I1 = setdiff(1:na, I0)
      if(length(I0)) l0[s] = max(loglik[I0])
      if(length(I1)) l1[s] = max(loglik[I1])
    }
  }
  GLR = exp(l0 - l1)
  data.frame(H0 = hyp, "loglik numerator" = l0, "loglik denominator"= l1,
             GLR = GLR)
}



