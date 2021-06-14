# binary type data

epsilon.binary<-2.0
check.binary.row<-function(arr)
{
  if (sum(!is.na(arr))==0)
  {
    return (F)
  }
  else
  {
    idx<-!is.na(arr)
    if (sum(arr[idx])==0 || sum(arr[idx])==sum(idx))
    {
      return (F)
    }
    else
    {
      return (T)
    }
  }
}
check.binary<-function(mat,name)
{
  index<-apply(mat,1,check.binary.row)
  n<-sum(!index)
  if (n>0)
  {
    w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
    warning(w)
  }
  mat_c<-mat[index,]
  rownames(mat_c)<-rownames(mat)[index]
  colnames(mat_c)<-colnames(mat)
  return (mat_c)
}

base.binary.row<-function(arr)
{
  idx<-!is.na(arr)
  n<-sum(idx)
  m<-sum(arr[idx])
  return (log(m/(n-m)))
}

base.binary<-function(mat)
{
  mat_b<-matrix(0,nrow(mat),ncol(mat))
  ar_b<-apply(mat,1,base.binary.row)
  mat_b[1:nrow(mat_b),]<-ar_b
  return (mat_b)
}

update.binary<-function(mat,mat_b,mat_now,eps)
{
  mat_p<-mat_b+mat_now
  mat_u<-matrix(0,nrow(mat),ncol(mat))
  idx1<-!is.na(mat) & mat==1
  idx0<-!is.na(mat) & mat==0
  index<-is.na(mat)
  arr<-exp(mat_p)
  mat_u[index]<-mat_now[index]
  mat_u[idx1]<-mat_now[idx1]+eps*epsilon.binary/(1.0+arr[idx1])
  mat_u[idx0]<-mat_now[idx0]-eps*epsilon.binary*arr[idx0]/(1.0+arr[idx0])
  return (mat_u)
}

stop.binary<-function(mat,mat_b,mat_now,mat_u)
{
  index<-!is.na(mat)
  mn<-mat_b+mat_now
  mu<-mat_b+mat_u
  arn<-exp(mn)
  aru<-exp(mu)
  idx1<-!is.na(mat) & mat==1
  idx0<-!is.na(mat) & mat==0
  lgn<-sum(log(arn[idx1]/(1+arn[idx1])))+sum(log(1/(1+arn[idx0])))
  lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
  return (lgu-lgn)
}

LL.binary<-function(mat,mat_b,mat_u)
{
  index<-!is.na(mat)
  mu<-mat_b+mat_u
  aru<-exp(mu)
  idx1<-!is.na(mat) & mat==1
  idx0<-!is.na(mat) & mat==0
  lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
  return (lgu)
}

LLmax.binary<-function(mat)
{
  return (0)
}

LLmin.binary<-function(mat,mat_b)
{
  index<-!is.na(mat)
  aru<-exp(mat_b)
  idx1<-!is.na(mat) & mat==1
  idx0<-!is.na(mat) & mat==0
  lgu<-sum(log(aru[idx1]/(1+aru[idx1])))+sum(log(1/(1+aru[idx0])))
  return (lgu)
}

binary_type_base <- function( data,dimension=2 ,name="test")
{
  data<-check.binary(data,name)
  data_b<-base.binary(data)
  data_now<-matrix(0,nrow(data),ncol(data))
  data_u<-update.binary(data,data_b,data_now)
  data_u<-nuclear_approximation(data_u,dimension)
  while (T)
  {
    thr<-stop.binary(data,data_b,data_now,data_u)
    message(thr)
    if (thr<0.2)
    {
      break
    }
    data_now<-data_u
    data_u<-update.binary(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
  }
  return (data_now)
}
