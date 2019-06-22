#library(mixer)
#library(blockmodels)

sbm_blockmodels<-function(mat,dist=c("Gauss","Bern","Poiss"),nmax,ef=1.05){
  # TODO : Check package for other distributions as well as covariates inclusion 
}

sbm_mixer<-function(mat,meth=c("classification","variational","bayesian"),Qmin=1,Qmax=50,Nbiter=50, Fpnbiter=25){
  icl=vector(mode='numeric',length=n)
  index=1:n
  mixture=mixer(mat,qmin=Qmin,qmax=Qmax,method = meth, nbiter=Nbiter, fpnbiter=Fpnbiter, improve=TRUE)
  return(mixture)
}

mixture_analysis<-function(mixture,n,criteria=c("max","relaxed","fix"),eps=0.05,fix=50,returnicl=TRUE){  ###n=(Qmax_Qmin)+1
  ### Find optimal number of groups according to criteria
  if(criteria=="fix"){
    Qoptim=min(n,fix)
  }
  else{
    icl=vector(mode='numeric',length=n)
    index=1:n
    
    for (i in index){
      ### get ICL per number of groups
      icl[i]=mixture$output[i][1][[1]]$criterion
    }
    
    Qmax=which.max(icl)
    iclmax=max(icl)
    
    if(criteria=="relaxed"){
      i=Qmax
      target=(1+eps)*iclmax
      
      while(icl[i]>target & i<=n){  ##> because it's negative
        i=i+1}
      Qoptim=i
    }
    else{
      Qoptim=Qmax}
  }
  
  ### Extract optimal label paritioning and group adjacency matrix
  adjacency=data.matrix(mixture$output[Qoptim][1][[1]]$Pis)
  taus=data.matrix(mixbayes$output[i][1][[1]]$Taus)
  
  partition=data.frame(row.names = mixbayes$nnames[[1]])
  lv=split(taus, rep(1:ncol(taus), each = nrow(taus)))
  partition$labels=unlist(lapply(lv,FUN = function(x) which.max(x)))
  
  if(returnicl){
    result=list(part=partition,adj=adjacency,icl=icl) 
  }
  else{
    result=list(part=partition,adj=adjacency)
  }
  
  return(result)  ## List of partition array, adjacency matrix and eventually icl table to plot or cross-check
}


####################################### Unit tests ##########################################################
mixture_analysis_test<-function(){ ## To execute row by row
  setwd(".")
  test_data<-load("Unit_tests/data.RData")
  anal<-NULL
  anal<-mixture_analysis(mixbayes,n,"relaxed",0.05,100)
  corrplot::corrplot(anal$adj)
  plot(anal$icl)
}

