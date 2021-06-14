check.matrix.element<-function(x)
{
  if (!is.matrix(x))
  {
    return (T)
  }
  else
  {
    return (F)
  }
}

ncol.element<-function(x)
{
  return (ncol(x))
}

nrow.element<-function(x)
{
  return (nrow(x))
}

check<-function(mat,type,name)
{
  if (type=="binary")
  {
    return (check.binary(mat,name))
  }
  else if (type=="gaussian")
  {
    return (check.gaussian(mat,name))
  }
  else if (type=="poisson")
  {
    return (check.poisson(mat,name))
  }
  else
  {
    e<-paste("unknown type ",type,sep="")
    stop(e)
  }
}

LRAcluster<-function(data,types,dimension=2,names=as.character(1:length(data)))
{
  eps<-0.0
  if (!is.list(data))
  {
    stop("the input data must be a list!")
  }
  c<-sapply(data,check.matrix.element)
  if (sum(c)>0)
  {
    stop("each element of input list must be a matrix!")
  }
  c<-sapply(data,ncol.element)
  if (length(levels(factor(c)))>1)
  {
    stop("each element of input list must have the same column number!")
  }
  if (length(data)!=length(types))
  {
    stop("data and types must be the same length!")
  }
  nSample<-c[1]
  loglmin<-0
  loglmax<-0
  loglu<-0.0
  nData<-length(data)
  for (i in 1:nData)
  {
    data[[i]]<-check(data[[i]],types[[i]],names[[i]])
  }
  nGeneArr<-sapply(data,nrow.element)
  nGene<-sum(nGeneArr)
  indexData<-list()
  k=1
  for(i in 1:nData)
  {
    indexData[[i]]<- (k):(k+nGeneArr[i]-1)
    k<-k+nGeneArr[i]
  }
  base<-matrix(0,nGene,nSample)
  now<-matrix(0,nGene,nSample)
  update<-matrix(0,nGene,nSample)
  thr<-array(0,nData)
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      base[indexData[[i]],]<-base.binary(data[[i]])
      loglmin<-loglmin+LLmin.binary(data[[i]],base[indexData[[i]],]) 
      loglmax<-loglmax+LLmax.binary(data[[i]])
    }
    else if (types[[i]]=="gaussian")
    {
      base[indexData[[i]],]<-base.gaussian(data[[i]])
      loglmin<-loglmin+LLmin.gaussian(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.gaussian(data[[i]])
    }
    else if (types[[i]]=="poisson")
    {
      base[indexData[[i]],]<-base.poisson(data[[i]])
      loglmin<-loglmin+LLmin.poisson(data[[i]],base[indexData[[i]],])
      loglmax<-loglmax+LLmax.poisson(data[[i]])
    }
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="gaussian")
    {
      update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
    else if (types[[i]]=="poisson")
    {
      update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
    }
  }
  update<-nuclear_approximation(update,dimension)
  nIter<-0
  thres<-array(Inf,3)
  epsN<-array(Inf,2)
  while(T)
  {
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        thr[i]<-stop.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="gaussian")
      {
        thr[i]<-stop.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
      else if (types[[i]]=="poisson")
      {
        thr[i]<-stop.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],update[indexData[[i]],])
      }
    }
    nIter<-nIter+1
    thres[1]<-thres[2]
    thres[2]<-thres[3]
    thres[3]<-sum(thr)
    epsN[1]<-epsN[2]
    epsN[2]<-eps
    if (nIter>5)
    {
      if (runif(1)<thres[1]*thres[3]/(thres[2]*thres[2]+thres[1]*thres[3]))
      {
        eps<-epsN[1]+0.05*runif(1)-0.025
      }
      else
      {
        eps<-epsN[2]+0.05*runif(1)-0.025
      }
      if (eps< -0.7)
      {
        eps<- 0
        epsN<-c(0,0)
      }
      if (eps > 1.4)
      {
        eps<-0
          epsN<-c(0,0)
      }
    }
    if (sum(thr)<nData*0.2)
    {
      break
    }
    now<-update
    for (i in 1:nData)
    {
      if (types[[i]]=="binary")
      {
        update[indexData[[i]],]<-update.binary(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="gaussian")
      {
        update[indexData[[i]],]<-update.gaussian(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
      else if (types[[i]]=="poisson")
      {
        update[indexData[[i]],]<-update.poisson(data[[i]],base[indexData[[i]],],now[indexData[[i]],],exp(eps))
      }
    }
    update<-nuclear_approximation(update,dimension)
  }
  for (i in 1:nData)
  {
    if (types[[i]]=="binary")
    {
      loglu<-loglu+LL.binary(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="gaussian")
    {
      loglu<-loglu+LL.gaussian(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
    else if (types[[i]]=="poisson")
    {
      loglu<-loglu+LL.poisson(data[[i]],base[indexData[[i]],],update[indexData[[i]],])
    }
  }
  sv<-svd(update,nu=0,nv=dimension)
  coordinate<-diag(c(sv$d[1:dimension],0))[1:dimension,1:dimension]%*%t(sv$v)
  colnames(coordinate)<-colnames(data[[1]])
  rownames(coordinate)<-paste("PC ",as.character(1:dimension),sep="")
  ratio<-(loglu-loglmin)/(loglmax-loglmin)
  return (list("coordinate"=coordinate,"potential"=ratio))
}
