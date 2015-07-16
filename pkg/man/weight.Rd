\name{weight}
\alias{weight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model averaging weights
}
\description{
This function gives the model averaging weights of the Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic models.
}
\usage{
weight(inc,time,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions.
}
}

\value{
A vector with the model averaging weights 
}
\references{
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information-theoretic Approach,
2nd Edition, Springer-Verlag, New York, 2002.\cr
G. Claeskens, N. L. Hjort, Model selection and model averaging, Cambridge University Press, 2008.

}

\author{
Carlos Sebrango, Lizet Sanchez, Ziv Shkedy
}

\examples{
require(nlme)
require(NRAIA)
require(drc)

##Incidence Data example 1
inc1 <-c(1, 0, 2, 5, 2, 2, 3, 5, 7, 14, 17, 18, 20, 38, 50, 45, 57,
44, 12, 11, 10, 2, 3, 0, 3, 2, 2, 0, 2, 0, 0, 0, 1, 2, 2, 1)

##Incidence Data example 2
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 25, 
10, 5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)

## time vectors
tii<- c(1:36)
ti2<- c(1:34)

## Model averaging Weights for Incidence data example 1

weight(inc1,time=tii,start=list(list(alpha=375,k=2.38,
gamma=0.76,eta=16),list(alpha=375,gamma=1,eta=16),
list(alpha=380,eta=13,beta=7,n=5),list(alpha=380,eta=20,
beta=0,gamma=1),list(alpha=410,eta=12,beta=11,k=3),
list(alpha=380,beta=4,g=1,eta=13,k=15)))

##or
weight(inc1,time=tii)

## Model averaging Weights for Incidence data example 2

weight(inc2,ti2,start=list(list(alpha=355,k=1,gamma=1,
eta=14),list(alpha=355,gamma=1,eta=14), list(alpha=355,
eta=13,beta=10,n=6), list(alpha=355,eta=11,beta=20,
gamma=1),list(alpha=355,eta=12,beta=22,k=3), 
list(alpha=355,beta=15,g=1,eta=13,k=10)))

##or
weight(inc2,ti2)
}
