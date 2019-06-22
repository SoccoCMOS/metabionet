
infer_net<-function(mbmat,rules){
  adjmat=matrix(data=rbinom(10 * 10, 1, 0.5),nrow=10,ncol=10)
  return(adjmat)
}