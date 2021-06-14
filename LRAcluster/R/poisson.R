# poisson type data

epsilon.poisson<-0.5

check.poisson.row<-function(arr)
{
  if (sum(!is.na(arr))==0)
  {
    return (F)
  }
  else
  {
    idx<-!is.na(arr)
    if (sum(arr[idx]<0)>0)
    {
      return (F)
    }
    else
    {
      return (T)
    }
  }
}

check.poisson<-function(mat,name)
{
  w<-paste(name," is poisson type. Add 1 to all counts",sep="")
  message(w)
  index<-apply(mat,1,check.poisson.row)
  n<-sum(!index)
  if (n>0)
  {
    w<-paste("Warning: ",name," have ",as.character(n)," invalid lines!",sep="")
    warning(w)
  }
  mat_c<-mat[index,]+1
  rownames(mat_c)<-rownames(mat)[index]
  colnames(mat_c)<-colnames(mat)
  return (mat_c)
}

base.poisson.row<-function(arr)
{
  idx<-!is.na(arr)
  m<-sum(log(arr[idx]))
  n<-sum(idx)
  return(m/n)
}

base.poisson<-function(mat)
{
  mat_b<-matrix(0,nrow(mat),ncol(mat))
  ar_b<-apply(mat,1,base.poisson.row)
  mat_b[1:nrow(mat_b),]<-ar_b
  return (mat_b)
}

update.poisson<-function(mat,mat_b,mat_now,eps)
{
  mat_p<-mat_b+mat_now
  mat_u<-matrix(0,nrow(mat),ncol(mat))
  index<-!is.na(mat)
  mat_u[index]<-mat_now[index]+eps*epsilon.poisson*(log(mat[index])-mat_p[index])
  index<-is.na(mat)
  mat_u[index]<-mat_now[index]
  return (mat_u)
}

stop.poisson<-function(mat,mat_b,mat_now,mat_u)
{
  index<-!is.na(mat)
  mn<-mat_b+mat_now
  mu<-mat_b+mat_u
  lgn<-sum(mat[index]*mn[index]-exp(mn[index]))
  lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
  return (lgu-lgn)
}

LL.poisson<-function(mat,mat_b,mat_u)
{
  index<-!is.na(mat)
  mu<-mat_b+mat_u
  lgu<-sum(mat[index]*mu[index]-exp(mu[index]))
  return (lgu)
}

LLmax.poisson<-function(mat)
{
  index<-!is.na(mat)
  lgu<-sum(mat[index]*log(mat[index])-mat[index])
  return (lgu)
}

LLmin.poisson<-function(mat,mat_b)
{
  index<-!is.na(mat)
  lgu<-sum(mat[index]*mat_b[index]-exp(mat_b[index]))
  return (lgu)
}

poisson_type_base<-function(data,dimension=2,name="test")
{
  data<-check.poisson(data,name)
  data_b<-base.poisson(data)
  data_now<-matrix(0,nrow(data),ncol(data))
  data_u<-update.poisson(data,data_b,data_now)
  data_u<-nuclear_approximation(data_u,dimension)
  while(T)
  {
    thr<-stop.poisson(data,data_b,data_now,data_u)
    message(thr)
    if (thr<0.2)
    {
      break
    }
    data_now<-data_u
    data_u<-update.poisson(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
  }
  return (data_now)
}
