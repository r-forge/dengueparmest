\name{allmodel.m}
\alias{allmodel.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak.
}
\description{
This function gives the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. Also this function, when all the built-in models are used,  gives the AIC of each model, the model averaged weights and predicted incidence and cumulative cases.this function is used retrospectively, that is, when the outbreak occurred.
}
\usage{
## Object of the S3 class dengue
allmodel.m(inc,time,start=NULL,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model="all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters. By default, the initial values are provided by self-starting functions.
}
\item{model}{
The nonlinear model to be used for fitting the data. Built-in models are  "Richards"", "Logistic3P", "SigmEmax", "Gompertz", "Weibull" and "Logistic5P". If model="all" the parameter estimate will be done taking into account all built-in models via model averaging.
}
}

\value{
 An object with the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{AIC}{The AIC for each built-in model and model averaged}
\item{tTable}{A table with parameter estimates and t test. It is not available when all the model are used.}
\item{Weights}{Model averaged weights. It is not availabe when is used only one built-in model}
\item{FinalSize}{95\% confidence interval and point estimate of the final size of outbreak for each built-in model and model averaged estimate}
\item{TurningPoint}{95\% confidence interval and point estimate of the turning point of outbreak for each built-in model and model averaged estimate}
\item{Predict}{Predicted cumulative cases for each built-in model}
\item{PredictMA}{Predicted cumulative cases for model averaged}
\item{PredInc}{Predicted incidence for each built-in model}
\item{PredMAInc}{Predicted incidence for model averaged}
\item{function.type}{Name of the function}
\item{model.type}{models used to estimate}
Generic functions such as plot and summary have methods to show the results of the fit
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
inc1 <-c(1, 0, 2, 5, 2, 2, 3, 5, 7, 14, 17, 18, 20, 38, 50, 45,
57, 44, 12, 11, 10, 2, 3, 0, 3, 2, 2, 0, 2, 0, 0, 0, 1, 2, 2, 1)

##Incidence Data example 2
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 
25, 10, 5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)

## time vectors
tii<- c(1:36)
ti2<- c(1:34)
## Parameter estimate for Richards model (for Incidence data example 1)

allmodel.m(inc1,time=tii, start=list(alpha=375,k=2.38,gamma=0.76,eta=16),model="Richards")

## or
allmodel.m(inc1,time=tii, model="Richards")

## Parameter estimate for each built-in model and model averaged 
##estimate for final size and turning point of outbreak 
##(for Incidence data example 1)

allmodel.m(inc1,time=tii,start=list(list(alpha=375,k=2.38,gamma=0.76,eta=16),
list(alpha=375,gamma=1,eta=16),list(alpha=380,eta=13,beta=7,n=5),
list(alpha=380,eta=20,beta=0,gamma=1),list(alpha=410,eta=12,beta=11,k=3),
list(alpha=380,beta=4,g=1,eta=13,k=15)),model="all")

## or
allmodel.m(inc1,time=tii,model="all")

# Parameter estimate for 3P Logistic model (for Incidence data example 2)

allmodel.m(inc2,ti2, start=list(alpha=375,gamma=1,eta=16), model="logistic3P")

## or
allmodel.m(inc2,ti2,model="logistic3P")

## Parameter estimate for each built-in model and model averaged estimate 
##for final size and turning point of outbreak 
##for Incidence data example 2
allmodel.m(inc2,ti2,start=list(list(alpha=355,k=1,gamma=1,eta=14),
list(alpha=355,gamma=1,eta=14), list(alpha=355,eta=13,beta=10,n=6),
list(alpha=355,eta=11,beta=20,gamma=1),list(alpha=355,eta=12,beta=22,k=3),
list(alpha=355,beta=15,g=1,eta=13,k=10)),model="all")

## or
allmodel.m(inc2,ti2,model="all")

}
