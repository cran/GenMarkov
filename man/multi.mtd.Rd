\name{multi.mtd}
\alias{multi.mtd}
\title{
Estimation of Multivariate Markov Chains
}
\description{
This function estimates the Mixture Distribution Model (Raftery (1985)) for Multivariare Markov Chains. It considers Berchtold (2001) optimization algorithm for the parameters and estimates the probabilities transition matrices as proposed in Ching (2002).
}
\usage{
multi.mtd(y, deltaStop = 0.0001, is_constrained = TRUE, delta = 0.1)
}
\arguments{
  \item{y}{Matrix of categorical data sequences.
}
  \item{deltaStop}{the delta below which the optimization phases of the parameters stop.
}
  \item{is_constrained}{flag indicating wether the function will consider the usual set of constraints (usual set: \emph{TRUE}, new set of constraints: \emph{FALSE}).
}
  \item{delta}{the amount of change to increase/decrease in the parameters for each iteration of the optimization algorithm.
}
}
\value{
The function returns a list with the parameter estimates,  standard-errors, z-statistics, p-values and the value of the log-likelihood function, for each equation.
}
\references{
Raftery, A. E. (1985). A Model for High-Order Markov Chains. Journal of the Royal Statistical Society. Series B (Methodological), 47(3), 528-539. \url{http://www.jstor.org/stable/2345788}

Berchtold,  A.  (2001).  Estimation  in  the  Mixture  Transition  Distribution  Model. Journal of Time Series Analysis, 22(4), 379-397.\doi{10.1111/1467-9892.00231}

Ching, W. K., E. S. Fung, and M. K. Ng (2002). A multivariate Markov chain model for categorical data sequences and its applications in demand predictions. IMA Journal of Management Mathematics, 13(3), 187-199. \doi{10.1093/imaman/13.3.187}
}
\author{
Carolina Vasconcelos and Bruno Damasio
}
\note{See details of the optimization procedure in Berchtold (2001) <\doi{10.1111/1467-9892.00231}>.
}
\examples{
set.seed(1234)
s1 <- sample(c(1,2), 500, replace=TRUE)
s2 <- sample(c(1,2), 500, replace=TRUE)
multi.mtd(y = cbind(s1,s2))
}
