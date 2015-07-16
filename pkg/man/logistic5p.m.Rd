\name{logistic5p.m}
\alias{logistic5p.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
5P Logistic model fit
}
\description{
This function fits the five parameter logistic model (5P Logistic ):\cr\cr
C(t)=beta + (alpha-beta)/[1+(2^(1/g)-1).(t/eta)^-k]^g \cr\cr
Where \cr
C(t): Cumulative number of reported cases at time t.\cr
alpha:  Final outbreak size or maximum  case number or carrying capacity.\cr
eta: Turning point of the epidemic.\cr
g: Asymmetry factor. \cr
k: Shape factor. \cr
beta: Lower asymptote.
}
\usage{
logistic5p.m (inc,time,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values for the parameters in the 5P logistic model. The list must be in the following way, list(alpha=alpha0,beta=beta0,g=g0,eta=eta0,k=k0), where alpha0,  gamma0 and eta0 are the starting values of 5P logistic model parameters. By default, the initial values are provided by self-starting function.
}
}

\value{
An object of a class gnls representing the nonlinear model fit
}
\references{
J. Liao, R. Liu, Re-parameterization of five-parameter logistic function, Journal of Chemometrics 23 (5), pp 248-253, 2009.
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

## 5P logistic model fit for Incidence data example 1
logistic5p.m(inc1,time=tii,start=list(alpha=380,beta=4,g=1,eta=13,k=15))

##or
logistic5p.m(inc1,time=tii)

## 5P logistic model fit for Incidence data example 2
logistic5p.m(inc2,ti2,list(alpha=355,beta=15,g=1,eta=13,k=10))

##or
logistic5p.m(inc2,ti2)


}
