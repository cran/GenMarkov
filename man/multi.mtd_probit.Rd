\name{multi.mtd_probit}
\alias{multi.mtd_probit}
\title{
Mixture Transition Distribution - Probit Model
}
\description{
Estimation of Multivariate Markov Chains through the proposed model by Nicolau (2014). This model presents two attractive features: it is completely free of constraints, thereby facilitating the estimation procedure, and it is more precise at estimating the transition probabilities of a multivariate or higher-order Markov chain than the Raftery's MTD model.
}
\usage{
multi.mtd_probit(y, initial, nummethod='bfgs')
}
\arguments{
  \item{y}{Matrix of categorical data sequences
}
  \item{initial}{Vector of initial values
}
  \item{nummethod}{Numerical maximisation method, currently either "NR" (for Newton-Raphson), "BFGS" (for Broyden-Fletcher-Goldfarb-Shanno), "BFGSR"     (for the BFGS algorithm implemented in R), "BHHH" (for Berndt-Hall-Hall-Hausman), "SANN" (for Simulated ANNealing), "CG" (for Conjugate     Gradients), or "NM" (for Nelder-Mead). Lower-case letters (such as "nr" for Newton-Raphson) are allowed. The default method is "BFGS". For more details see \code{\link{maxLik}}.
}
}
\value{
The function returns a list with the parameter estimates,  standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
}
\references{
Nicolau, J. (2014). A new model for multivariate markov chains. Scandinavian Journal of Statistics, 41(4), 1124-1135.\doi{10.1111/sjos.12087}
}
\author{
Carolina Vasconcelos and Bruno Damasio
}
\examples{
set.seed(1234)
s1 <- sample(c(1,2), 500, replace=TRUE)
s2 <- sample(c(1,2), 500, replace=TRUE)
multi.mtd_probit(y = cbind(s1,s2), initial=c(1,1,1), nummethod='bfgs')
}

