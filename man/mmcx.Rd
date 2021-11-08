\name{mmcx}
\alias{mmcx}

\title{
Non-homogeneous Multivariate Markov Chains
}
\description{
Estimates Multivariate Markov Chains that depend on a exogeneous variables. The model is based on the Mixture Transition Distribution model, and considers non-homogeneous Markov Chains, instead of homogeneous Markov Chains as in Raftery (1985).
}
\usage{
mmcx(y,x,initial)
}
\arguments{
  \item{y}{Matrix of categorical data sequences.
}
  \item{x}{Matrix of covariates (exogeneous variables).
}
  \item{initial}{Vector of initial values.
}
}
\value{
The function returns a list with the parameter estimates,  standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
}
\references{
Raftery, A. E. (1985). A Model for High-Order Markov Chains. Journal of the Royal Statistical Society. Series B (Methodological), 47(3), 528-539. \url{http://www.jstor.org/stable/2345788}

Ching, W. K., E. S. Fung, and M. K. Ng (2002). A multivariate Markov chain modelfor categorical data sequences and its applications in demand predictions. IMA Journal of Management Mathematics, 13(3), 187-199. \doi{10.1093/imaman/13.3.187}
}
\author{
Carolina Vasconcelos and Bruno Damasio
}
\seealso{
Optmization is done through \code{\link{auglag}}.
}
\examples{
set.seed(1234)
s1 <- sample(c(1,2), 500, replace=TRUE)
s2 <- sample(c(1,2), 500, replace=TRUE)
x <- rnorm(500)
mmcx(y = cbind(s1,s2), x = cbind(x), initial=c(1,1))
}

