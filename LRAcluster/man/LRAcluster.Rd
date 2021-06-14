\name{LRAcluster}
\alias{LRAcluster}
\title{integrated analysis of cancer omics data by low-rank approximation}
\description{
The LRAcluster function is the main interface of this package, it gets a list of matrix as input and outputs the coordinates of the samples in the reduced space and the explained potential.
}
\usage{
LRAcluster(data,types,dimension=2,names)
}
\arguments{
  \item{data}{a list of data matrix,please keep the columns are the same order of samples}
  \item{types}{a list of data types, can be binary, gaussian, or poisson}
  \item{dimension}{the reduced dimension}
  \item{names}{data names}
}
\seealso{
  \code{\link{testdata}}
}
\value{
  \item{coordinate}{A matrix of the coordinates of all the samples in the reduced space}
  \item{potential}{ratio of explained variance}
}
\examples{
data(testdata)
cols<-array(NA,length(labels))
cols[labels=="BRCA"]<-"red"
cols[labels=="COAD"]<-"green"
cols[labels=="LUAD"]<-"blue"
rlist<-LRAcluster(data=data,types=types,names=names,dimension=2)
rclust<-kmeans(t(rlist$coordinate),3)
plot(rlist$coordinate[1,],rlist$coordinate[2,],xlab="pc 1",ylab="pc 2",col=cols,main="demo")
legend("topright",legend=c("BRCA","COAD","LUAD"),pch=1,col=c("red","green","blue"))
}
