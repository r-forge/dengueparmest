\name{plot}
\alias{plot.dengue}
\title{
Plots for a visualization of the results
}
\description{
plot is a generic function used to produce plots for a visualization of the results  of the functions: allmodel.m, allmodel.predict and changetimeFSTP from the S3 class dengue. The function invokes particular methods of the S3 class dengue. 
}
\usage{
\method{plot}{dengue}(x, which = c(1, 2), xlab = "", \dots)
}
\arguments{
  \item{x}{
An object for which a plot is desired. It could be any of the following objects: allmodel.m, allmodel.predict and changetimeFSTP
}
  \item{which}{
Different graph options, from 1 to 6. 1 draws a graph with cumulative number of cases and fitted models for each built-in models and model averaged, 2 draws a graph with incidence and fitted models, 3 draws a plot with the final size 95\% confidence interval and point estimate for each
built-in model and model averaged estimate, 4 draws a plot with the turning point 95\% confidence interval and point estimate for each built-in model and model averaged estimate, 5 draws a plot with the change over time of the final size estimates for each built-in models and model averaged and option 6 draws a plot with the change over time of the turning point estimates for each built-in models and model averaged.
}
  \item{xlab}{
A title for the x axis. It is not necesarry for options 3, 4, 5 y 6.
}
  \item{\dots}{
Additional arguments affectingthe plot produced. Arguments to be passed to methods, such as graphical parameters(legend, lines, etc) 
}
}
\details{
When all the built-in models and model averaged are used in the options 1, 2, 5 and 6 the plot displays all fitted models with the following colours: (Richards (red),3P logistic (green), Sigmoid Emax (dark blue), Gompertz (clear blue),Weibull (pink) and 5P logistic (yelow))and the fitted model via model averaging (gray colour). \cr
When the objects are allmodel.m and allmodel.predict can be used only the options 1, 2, 3 and 4. When the object is changetimeFSTP can only be used the options 5 and 6.
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
require(ggplot2)
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

## plot function for a allmodel.m object 
##Richards model (Incidence data example 1)
ee<-allmodels(inc1,tii, model="Richards")
plot(ee,which=c(1,2))

##using all built-in models (Incidence data example 2)
ee2<-allmodels(inc2,ti2, model="all")
plot(ee2,which=c(1,2,3,4))

##plot function for allmodel.predict object
##3P logistic model (Incidence data example 1)
## information avalaible to the time point 20
## observation to predict: 30 
ee3<-allmodelpredict(inc2[1:20],ti2[1:20],30, model="logistic3P")
plot(ee3,which=c(1,2))

##using all built-in models (Incidence data example 1)
## information avalaible to the time point 22
## observation to predict: 30 
ee4<-allmodelpredict(inc2[1:22],ti2[1:22],30, model="all")
plot(ee4,which=c(1,2,3,4))

##plot function for changetimeFSTP object
##Time point from which is required to see the changes over time 
##of the final size and turning point estimates:25
##(Incidence data example 1)
t1<-changetimeFSTP(inc1,tii,25)
plot(t1,which=c(5,6))

##Time point from which is required to see the changes over time 
##of the final size and turning point estimates:20
##(Incidence data example 2)
t2<-changetimeFSTP(inc2,ti2,20)
plot(t2,which=c(5,6))

}

