\name{summary}
\alias{summary.dengue}
\title{
Result summaries of the functions
}
\description{
summary is a generic function used to produce result summaries of the functions: allmodel.m, allmodel.predict and changetimeFSTP from the S3 class dengue. The function invokes particular methods of the S3 class dengue. 
}
\usage{
\method{summary}{dengue}(object, \dots)
}

\arguments{
  \item{object}{
An object for which a summary is desired. It could be any of the following objects: allmodel.m, allmodel.predict and changetimeFSTP.
}
\item{\dots}{
Additional arguments affecting the summary produced.
}
}

\value{
The form of the value returned by summary depends on the object of the S3 class dengue. In the case of allmodel.m returns the AIC, model averaged Weights and model averaged estimates for the final size and turning point of outbreak. If the object is allmodel.predict to the returned value is added a prediction of the incidence and cumulative cases in the required future observation. When the object is changetimeFSTP function the returned values are two tables with the changes over time of the final size and turning point for the time period required.
}
\references{
K. Burnham, D. R. Anderson, Model Selection and Multimodel Inference: A Practical Information-theoretic Approach,
2nd Edition, Springer-Verlag, New York, 2002.\cr
G. Claeskens, N. L. Hjort, Model selection and model averaging, Cambridge University Press, 2008.

}
\author{
Carlos Sebrango, Lizet Sanchez, Ziv Shkedy, Ewoud DE Troyer
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

## summary function for a allmodel.m object 
##Richards model (Incidence data example 1)
ee<-allmodels(inc1,tii, model="Richards")
summary(ee)

##using all built-in models (Incidence data example 2)
ee2<-allmodels(inc2,ti2, model="all")
summary(ee2)

##summary function for allmodel.predict object
##3P logistic model (Incidence data example 1)
## information avalaible to the time point 20
## observation to predict: 30 
ee3<-allmodelpredict(inc2[1:20],tii[1:20],30, model="logistic3P")
summary(ee3)

##using all built-in models (Incidence data example 1)
## information avalaible to the time point 22
## observation to predict: 30 
ee4<-allmodelpredict(inc2[1:22],tii[1:22],30, model="all")
summary(ee4)

##summary function for changetimeFSTP object
##Time point from which is required to see the changes over time 
##of the final size and turning point estimates:25
##(Incidence data example 1)
t1<-changetimeFSTP(inc1,tii,25)
summary(t1)

##Time point from which is required to see the changes over time 
##of the final size and turning point estimates:20
##(Incidence data example 2)
t2<-changetimeFSTP(inc2,ti2,20)
summary(t2)


}


