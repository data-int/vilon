# gaussian type data
epsilon.gaussian=0.5

check.gaussian.row<-function(arr)
{
  if (sum(!is.na(arr))==0)
  {
    return (F)
  }
  else
  {
    return (T)
  }
}
check.gaussian<-function(mat,name)
{
  index<-array(T,nrow(mat))
  for(i in 1:nrow(mat))
  {
    if (sum(is.na(mat[i,])==ncol(mat)))
    {
      war<-paste("Warning: ",name,"'s ",as.character(i)," line is all NA. Delete this line",sep="")
      warning(war)
      index[i]<-F
    }
  }
  mat_c<-mat[index,]
  rownames(mat_c)<-rownames(mat)[index]
  colnames(mat_c)<-colnames(mat)
  return (mat_c)
}

base.gaussian.row<-function(arr)
{
  idx<-!is.na(arr)
  return (mean(arr[idx]))
}

base.gaussian<-function(mat)
{
  mat_b<-matrix(0,nrow(mat),ncol(mat))
  ar_b<-apply(mat,1,base.gaussian.row)
  mat_b[1:nrow(mat_b),]<-ar_b
  return (mat_b)
}

update.gaussian<-function(mat,mat_b,mat_now,eps)
{
  mat_p<-mat_b+mat_now
  mat_u<-matrix(0,nrow(mat),ncol(mat))
  index<-!is.na(mat)
  mat_u[index]<-mat_now[index]+eps*epsilon.gaussian*(mat[index]-mat_p[index])
  index<-is.na(mat)
  mat_u[index]<-mat_now[index]
  return (mat_u)
}

stop.gaussian<-function(mat,mat_b,mat_now,mat_u)
{
  index<-!is.na(mat)
  mn<-mat_b+mat_now
  mu<-mat_b+mat_u
  ren<-mat[index]-mn[index]
  reu<-mat[index]-mu[index]
  lgn<- -0.5*sum(ren*ren)
  lgu<- -0.5*sum(reu*reu)
  return (lgu-lgn)
}

LL.gaussian<-function(mat,mat_b,mat_u)
{
  index<-!is.na(mat)
  mu<-mat_b+mat_u
  reu<-mat[index]-mu[index]
  lgu<- -0.5*sum(reu*reu)
  return (lgu)
}

LLmax.gaussian<-function(mat)
{
  return (0.0)
}

LLmin.gaussian<-function(mat,mat_b)
{
  index<-!is.na(mat)
  reu<-mat[index]-mat_b[index]
  lgu<- -0.5*sum(reu*reu)
  return (lgu)
}

gaussian_base<-function(data,dimension=2,name="test")
{
  data<-check.gaussian(data,name)
  data_b<-base.gaussian(data)
  data_now<-matrix(0,nrow(data),ncol(data))
  data_u<-update.gaussian(data,data_b,data_now)
  data_u<-nuclear_approximation(data_u,dimension)
  while(T)
  {
    thr<-stop.gaussian(data,data_b,data_now,data_u)
    message(thr)
    if (thr<0.2)
    {
      break
    }
    data_now<-data_u
    data_u<-update.gaussian(data,data_b,data_now)
    data_u<-nuclear_approximation(data_u,dimension)
  }
  return (data_now)
}
