% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testLR.R
\name{testLR}
\alias{testLR}
\title{Finds classical likelihood ratio test statistics}
\usage{
testLR(pedigrees, sim, truePed = NULL, ref = 1, H0 = 1)
}
\arguments{
\item{pedigrees}{a list of pedigrees}

\item{sim}{pedigree with simulated marker data}

\item{truePed}{index or name of pedigree with marker data}

\item{ref}{index or name of pedigree in the denominator of each LR}

\item{H0}{index or name of pedigree consistent with H0}
}
\value{
A data frame with LRs, the max and log likelihood of the first pedigree,
the reference and optionally a plot.
}
\description{
n named pedigrees are given. For one of this we have simulated marker data.
The test statistic Delta = -ln(R) is calculated. Here is the ratio of the
maximum likelihood (or LR) for pedigrees consistent
with H0 to to the overall maximum.
1 to n,..., (n-1) to n and the maximum Z = max(LR.1.n, ..., LR.n-1.n)
are calculated
}
\examples{
pedigrees = examplePedigrees(plot = F, gm = F)
names(pedigrees) = paste0("Ped", 0:3)
truePed  = "Ped3"
ref = "Ped0"
nsim = 100
conditional = F
pedigrees[[truePed]] = setMarkers(pedigrees[[truePed]],
                  locusAttributes = NorwegianFrequencies[1:22])
pedigrees[[truePed]] = setMutmod(pedigrees[[truePed]], model = "proportional", rate = 0.001)
sim = simPedigree(pedigrees[[truePed]], idTarget ="2", idReferences = c("1", "3"),
                  nsim = nsim, seed = 1729, conditional = conditional)
H0 = c(1,4)
H0 = c("Ped0", "Ped3")
res = testLR(pedigrees, sim, truePed,  ref, H0)
res

}
