nuclear_approximation<-function(mat,dimension)
{
  svd<-svd(mat,nu=0,nv=0)
  if (dimension<length(svd$d))
  {
    lambda<-svd$d[dimension+1]
    svd<-svd(mat,nu=dimension,nv=dimension)
    indexh<-svd$d>lambda
    indexm<-svd$d<lambda
    dia<-array(svd$d,length(svd$d))
    dia[indexh]<-dia[indexh]-lambda
    dia[indexm]<-0
    mat_low<-svd$u%*%diag(c(dia[1:dimension],0))[1:dimension,1:dimension]%*%t(svd$v)
  }
  else
  {
    mat_low<-mat
  }
  return (mat_low)
}
