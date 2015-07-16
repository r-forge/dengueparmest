\name{weibull.m}
\alias{weibull.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weibull model fit
}
\description{
This function fits the Weibull model:\cr\cr
C(t)=beta + (alpha-beta).exp(-(t/eta)^-k)\cr\cr
Where \cr
C(t): Cumulative number of reported cases at time t.\cr
alpha:  Final outbreak size or maximum  case number or carrying capacity.\cr
eta: Turning point of the epidemic.\cr
k: Shape factor.\cr
beta: Lower asymptote.
}
\usage{
weibull.m (inc,time,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values for the parameters in the Weibull model. The list must be in the following way, list(alpha=alpha0,eta=eta0,beta=beta0,k=k0), where alpha0,eta0, beta0 and k0 are the starting values of Weibull  model parameters. By default, the initial values are provided by self-starting function.
}
}

\value{
An object of a class gnls representing the nonlinear model fit
}
\references{
G. Seber, C. Wild, Nonlinear regression, Wiley, New York, 1989.
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

## Weibull model fit for Incidence data example 1
weibull.m(inc1,time=tii,list(alpha=410,eta=12,beta=11,k=3))

##or
weibull.m(inc1,time=tii)

## Weibull model fit for Incidence data example 2
weibull.m(inc2,ti2,list(alpha=355,eta=12,beta=22,k=3))

##or
weibull.m(inc2,ti2)


}
