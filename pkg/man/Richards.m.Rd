\name{Richards.m}
\alias{Richards.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Richards model fit
}
\description{
This function fits the Richards model:\cr\cr
C(t)=alpha/[1+k.exp(-k.gamma.(t-eta))]^(1/k)\cr\cr
Where \cr
C(t): Cumulative number of reported cases at time t.\cr
alpha:  Final outbreak size or maximum  case number or carrying capacity.\cr
gamma:  Per capita growth rate or intrinsic growth rate.\cr
eta: Turning point of the epidemic.\cr
k: Exponent of the deviation from standard logistic curve.

}
\usage{
Richards.m(inc,time,start=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values for the parameters in the Richards model. The list must be in the following way, list(alpha=alpha0, k=k0, gamma=gamma0, eta=eta0), where alpha0, k0, gamma0 and eta0 are the starting values of Richards model parameters. By default, the initial values are provided by self-starting function.
}
}

\value{
An object of a class gnls representing the nonlinear model fit
}
\references{
F. Richards, A flexible growth function for empirical use, Journal of Experimental Botany 10 (29), pp 290-301, 1959.\cr
Y.H. Hsieh, Temporal trend and regional variability of 2001-2002 multiwave DENV-3 epidemic in Havana City: did Hurricane Michelle contribute to its severity?, Tropical Medicine and International Health, Vol. 18, no. 7, pp 830-838, 2013.\cr
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
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 25, 
10, 5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)

## time vectors
tii<- c(1:36)
ti2<- c(1:34)

## Richads model fit for Incidence data example 1
Richards.m(inc1,time=tii,start=list(alpha=375,k=2.38,gamma=0.76,eta=16))

##or
Richards.m(inc1,time=tii)

## Richads model fit for Incidence data example 2
Richards.m(inc2,ti2,list(alpha=355,k=1,gamma=1,eta=14))

##or
Richards.m(inc2,ti2)


}
