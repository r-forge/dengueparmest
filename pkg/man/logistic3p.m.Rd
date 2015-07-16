\name{logistic3p.m}
\alias{logistic3p.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3P Logistic model fit
}
\description{
This function fits the three parameter logistic model (3P Logistic ) proposed by Verhults:\cr\cr
C(t)=alpha/[1+exp(-gamma.(t-eta))]\cr\cr
Where \cr
C(t): Cumulative number of reported cases at time t.\cr
alpha:  Final outbreak size or maximum  case number or carrying capacity.\cr
gamma:  Per capita growth rate or intrinsic growth rate.\cr
eta: Turning point of the epidemic.
}
\usage{
logistic3p.m (inc,time,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values for the parameters in the 3P logistic model. The list must be in the following way, list(alpha=alpha0, gamma=gamma0, eta=eta0), where alpha0,  gamma0 and eta0 are the starting values of 3P logistic model parameters. By default, the initial values are provided by self-starting function.
}
}

\value{
An object of a class gnls representing the nonlinear model fit
}
\references{
A. Tsoularis, J. Wallace, Analysis of logistic growth models, Mathematical Biosciences, Vol. 179, no. 1, pp 21-55, 2002.
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
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 25, 10,
5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)

## time vectors
tii<- c(1:36)
ti2<- c(1:34)

## 3P logistic model fit for Incidence data example 1
logistic3p.m(inc1,time=tii,start=list(alpha=375,gamma=0.76,eta=16))

##or
logistic3p.m(inc1,time=tii)

## 3P logistic model fit for Incidence data example 2
logistic3p.m(inc2,ti2,list(alpha=355,gamma=1,eta=14))

##or
logistic3p.m(inc2,ti2)


}
