\name{changetimeFSTP}
\alias{changetimeFSTP}
\title{
Changes over time of the final size and turning point estimates for each built-in model and model averaged.
}
\description{
This function provides the changes over time of the final size and turning point estimates for each built-in model and model averaged. The output of the function are two tables, one for final size and the other for turning point. In the row appear the different time periods since is required until the end of time point available. In the columns are the built-in model and model averaged.
}
\usage{
## Object of the S3 class dengue
changetimeFSTP(inc, time, ini, start = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
 \item{ini}{
Time point from which is required to see the changes over time of the final size and turning point estimates
}
  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model="all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions.
}
}
\value{
 An object with the changes over time of the final size and turning point estimates for each built-in model and model averaged. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{Period}{Time period required to see the changes over time of the final size and turning point estimates}
\item{FSchangetime}{A Table with the changes over time of the final size for the time period required}
\item{TPchangetime}{A Table with the changes over time of the turning point for the time period required}
\item{function.type}{Name of the function}
\item{model.type}{models used to estimate}
Generic functions such as plot and summary have methods to show the results of the fit.
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
require(drc)
require(NRAIA)


##Incidence Data example 1
inc1 <-c(1, 0, 2, 5, 2, 2, 3, 5, 7, 14, 17, 18, 20, 38, 50, 45,
57, 44, 12, 11, 10, 2, 3, 0, 3, 2, 2, 0, 2, 0, 0, 0, 1, 2, 2, 1)

##Incidence Data example 2
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 
25, 10, 5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)

## time vectors
tii<- c(1:36)
ti2<- c(1:34)

## (for Incidence data example 1)
## Changes over time, since time point 25 to time point 36, of the final size and turning point 
## estimates for each built-in model and model averaged

changetimeFSTP(inc1,tii,25)


## (for Incidence data example 2)
## Changes over time, since time point 19 to time point 34, of the final size and turning point 
## estimates for each built-in model and model averaged

changetimeFSTP(inc2,ti2,19)
}

