\name{allmodel.predict}
\alias{allmodel.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Real-time predictions of the final size and the turning point at the end of the epidemic, the incidence and the cumulative number of cases in future observations.
}
\description{
This function provides real-time predictions of the final size and the turning point at the end of the epidemic for each built-in model and model averaged, as well as the incidence and the cumulative number of cases in future observations. Also this function, when all the built-in models are used,  gives the AIC of each model, the model averaged weights and predicted incidence and cumulative cases.
}
\usage{
## Object of the S3 class dengue
allmodel.predict(inc,time,pred,start=NULL,model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inc,time}{
Vector of equal length specifying incidence (number of reported cases per time unit) and time interval (from the start of outbreak).
}
\item{pred}{
Number of observation in which the incidence and the cumulative number of cases will be predicted.
}

  \item{start}{
A list with the starting values of the model to be used for fitting the data. If model="all" the imput must be a list of a list with the starting values of Richards, 3P logistic, Sigmoid Emax, Gompertz, Weibull and 5P logistic model parameters.By default, the initial values are provided by self-starting functions.
}
\item{model}{
The nonlinear model to be used for fitting the data. Built-in models are  "Richards"", "Logistic3P", "SigmEmax", "Gompertz", "Weibull" and "Logistic5P". If model="all" the parameter estimate will be done taking into account all built-in models via model averaging.
}
}

\value{
 An object with the parameter estimate for each built-in model and model averaged estimate for final size and turning point of outbreak. It is a list:
\item{Incidence }{All the available incidences}
\item{Time}{All the available time points}
\item{PredTime}{Period of time for which the prediction is required}
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


#Incidence Data example 1
inc1 <-c(1, 0, 2, 5, 2, 2, 3, 5, 7, 14, 17, 18, 20, 38, 50, 45,
57, 44, 12, 11, 10, 2, 3, 0, 3, 2, 2, 0, 2, 0, 0, 0, 1, 2, 2, 1)
#Incidence Data example 2
inc2<-c(5, 4, 4, 3, 7, 1, 4, 4, 9, 21,19, 33, 49, 45, 46, 28, 
25, 10, 5, 6, 2, 1, 4, 5, 0, 3, 2, 2, 1, 0, 3, 0, 2, 2)
# time vectors
tii<- c(1:36)
ti2<- c(1:34)

## Using only the information until time point 20 (data example 1)
# Prediction of the final size and turning point at the end of epidemic, the incidence and the cumulative number of cases in the observation number 36 using only the Richards model 

p1<-allmodel.predict(inc2[1:20],tii[1:20],36, model="Richards")
summary(p1)

## Using only the information until time point 22
## Now using all built-in model, Prediction of the final size and turning point at the end of ## epidemic, the incidence and the cumulative number of cases in the observation number 30 

p2<-allmodel.predict(inc2[1:22],tii[1:22],30, model="all")
summary(p2)

## (data example 2)
## Using only the information until time point 18, Prediction of the final size and turning 
## point at the end of epidemic, the incidence and the cumulative number of cases in the 
## observation number 31 using only the 3P logistic model

p3<-allmodel.predict(inc2[1:18],ti2[1:18],31, model="logistic3P")
summary(p3)

## Using only the information until time point 20
## Now using all built-in model, Prediction of the final size and turning point at the end of ## epidemic, the incidence and the cumulative number of cases in the observation number 30 

p4<-allmodel.predict(inc2[1:20],ti2[1:20],34, model="all")
summary(p4)

}