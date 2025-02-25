% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MoNET.R
\name{BuildNet}
\alias{BuildNet}
\title{Multi-omic network construction}
\usage{
BuildNet(threshold = 0)
}
\arguments{
\item{threshold}{A constant, gives a threshhold for removing nodes whose gene interaction score less than it, with the default value 0.}
}
\value{
a dataframe of network with two columns
}
\description{
Remove nodes whose gene interaction score less than the threshold from the muiti-omic network we manually integrated.
}
\examples{
MulOmicNet<-BuildNet(threshold=0.5)

}
