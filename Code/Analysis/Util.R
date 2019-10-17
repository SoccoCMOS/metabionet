############### Set of utility functions to help compute metrics on food webs #################
util_haslinkto<-function(net,targets,dir=c("in","out")){  ###returns lists of taxa labels that are connected to targets => example here targets are all groups of resources
  ###net: igraph object representing the network of interactions 
  ###targets: nodes from where to start energy pathways 
  ###direction of energy: in (they consume) or out (they are consumed)
  vtargets=V(net)[name %in% targets]
  l=list()
  for (v in targets){
    l=append(l,neighbors(net,v,dir)$name)
  }
  
  return(unlist(l))
}

containsAttr<-function(xset,yset){
  #print("New")
  #print(x)
  #print(y)
  return(length(intersect(xset,yset))>0)
}

get_energy_roots<-function(net,attribute=c("broad","scName"),target_lists,filter_attribute=c("broad","scName"),filter_values){
  ###net: igraph object representing the network of interactions
  ###attribute: node attribute from which we subset 
  ###target_lists: attribute values to subset
  #vertex_attr(net, attribute)%in%target_lists
  vroots=V(net)[
    unlist(lapply(vertex_attr(net, attribute),function(x) containsAttr(x,target_lists)))
    #vertex_attr(net, attribute)%in%target_lists
  ]
  
  haslinktovroots=unique(util_haslinkto(net,vroots,"out"))
  if(length(filter_values)>0){
    cond1=unlist(lapply(vertex_attr(net, filter_attribute),function(x) containsAttr(x,filter_values)))
    cond2=vertex_attr(net, "name")%in%haslinktovroots
    starters=V(net)[cond1 & cond2]
  }else{
    starters=V(net)[name %in% haslinktovroots]
  }
  return(starters)
}

pathways_stat<-function(net,node_sources){    
### Computes pathways statistics => 
  #net: igraph object
  ## node_sources: output of get_energy_roots 
  if(vcount(net)==0){
    print("Warning empty or disconnected network")
    brg_stat=list(sum_pathsize=0,mean_pathsize=0,median_pathsize=0,sd_pathsize=0,covered_taxa=NULL,nbcovered_taxa=0)
  }
  else{
    deb <<- net
    sums=vector(mode="numeric",length=length(node_sources))
    moys=vector(mode="numeric",length=length(node_sources))
    meds=vector(mode="numeric",length=length(node_sources))
    sds=vector(mode="numeric",length=length(node_sources))
    maxs=vector(mode="numeric",length=length(node_sources))
    
    ncov=vector(mode="character")
    
    cpt=0
    for (root in node_sources){
      cpt=cpt+1
      d=bfs(net, root, neimode = c("out"), 
            unreachable = FALSE, restricted = NULL, order = TRUE, rank = FALSE,
            father = FALSE, pred = FALSE, succ = FALSE, dist = TRUE,
            callback = NULL, extra = NULL, rho = parent.frame()) 
      
      roads=d$dist[which(!is.na(d$dist))]
      sums[cpt]=sum(roads)
      moys[cpt]=mean(roads)
      meds[cpt]=median(roads)
      sds[cpt]=sd(roads)
      maxs[cpt]=max(roads)
      ncov=union(ncov,names(which(!is.na(d$dist))))
    }
    
    nbcov=length(ncov)
    brg_stat=list(sum_pathsize=sum(sums),mean_pathsize=mean(moys),median_pathsize=median(meds),sd_pathsize=sd(sds),covered_taxa=ncov,nbcovered_taxa=nbcov,max_dist=max(maxs))
  }
  return(brg_stat)
}

ratio<-function(a,b){
  if(is.na(a) || is.na(b)){
    rat=0
  }
  else{
    if(b==0){
      rat=0
    }
    else{
      rat=a/b
    }
  }
}

energypathways_detail <- function(net,x,y,labx="x",laby="y"){   ### Analysis of a single network
  X=pathways_stat(net,x)
  Y=pathways_stat(net,y)
  
  ### X/Y statistics
  XF=data.frame(source=labx,sum=X$sum_pathsize, med=X$median_pathsize,mean=X$mean_pathsize,sd=X$sd_pathsize,nbcovered=X$nbcovered_taxa)
  YF=data.frame(source=laby,sum=Y$sum_pathsize, med=Y$median_pathsize,mean=Y$mean_pathsize,sd=Y$sd_pathsize,nbcovered=Y$nbcovered_taxa)
  comb=data.frame(source=paste(labx,laby,sep="/"),sum=ratio(X$sum_pathsize,Y$sum_pathsize), med=ratio(X$median_pathsize,Y$median_pathsize),mean=ratio(X$mean_pathsize,Y$mean_pathsize),sd=ratio(X$sd_pathsize,Y$sd_pathsize),nbcovered=ratio(X$nbcovered_taxa,Y$nbcovered_taxa))
  out=list(X=XF,Y=YF,RatioXY=comb)
  
  return(out)
}

energy_pathways_ut=function(){
  phototrophs=get_energy_roots(net,attribute="broad",target_lists=broad_light,filter_attribute="",filter_values=c())
  decomposers=get_energy_roots(net,attribute="broad",target_lists=broad_som,filter_attribute="",filter_values=c())
  brown_green_stats=energypathways_detail(net,y=phototrophs,x=decomposers,laby="green",labx="brown")
  
  bacterias=get_energy_roots(net,attribute="broad",target_lists=broad_som,filter_attribute="broad",filter_values=broad_bacteria)
  fungis=get_energy_roots(net,attribute="broad",target_lists=broad_som,filter_attribute="broad",filter_values=broad_fungi)
  bact_fungi_stats=energypathways_detail(net,y=bacterias,x=fungis,laby="bacterias",labx="fungis")
}

aggregate_occur=function(occur,cols,groups){
  othercols=setdiff(colnames(occur),cols)
  agg=data.frame(occur[,othercols])
  colnames(agg)=othercols
  new_names=unique(groups)
  for(g in new_names){
    c=names(which(groups==g))
    agg[,g]<-rowSums(occur[,c])
  }
  
  return(list(df=agg,names=new_names))
}

troph_sbm<-function(metaweb,file_raw_adjacency,file,file_sbm_groups,file_sbm_adjacency,minim=5,ef=1.5,maxim=Inf){
  ### Clustering ###
  #weighted_adj=data.frame(metaweb$get_adjacency_matrix("weight",F))
  bin_adj=data.frame(get.adjacency(metaweb$graph,
                                   attr=NULL, names=TRUE, sparse=FALSE))
  
  nodelist=data.frame(keys=unlist(V(metaweb$graph)$name))
  nodenames=data.frame(names=unlist(V(metaweb$graph)$scName))
  
  # Saving #
  #write.csv2(weighted_adj,file_raw_weighted_adjacency)
  write.csv2(bin_adj,file_raw_adjacency)
  
  sbm_model <- BM_bernoulli(
    membership_type="SBM",
    adj=as.matrix(bin_adj),
    verbosity=6,
    autosave='checkpoints_bin',
    plotting=character(0),
    exploration_factor=ef,
    explore_min=minim,
    explore_max=maxim,
    ncores=detectCores())
  
  sbm_model$estimate()
  argmax=which.max(sbm_model$ICL)
  
  ### Number of clusters detected by SBM ###
  cat("Number of clusters detected by SBM (maxICL) ",argmax) 
  
  ### Getting SBM parameters for best estimated model ###
  sbm_adj=sbm_model$model_parameters[[argmax]]$pi
  sbm_mem=sbm_model$memberships[[argmax]]
  sbm_mem$plot()
  trophic_groups=data.frame(taxa_key=nodelist,taxa_names=nodenames,trophic_group=apply(sbm_mem$Z,1,which.max))
  
  write.csv2(trophic_groups,file=file_sbm_groups)
  write.csv2(sbm_adj,file_sbm_adjacency)  
  
  return(list(tg=trophic_groups,adj=sbm_adj))
}

betanet<-function(gList,groups=NULL,pt_names=NULL,div='P',file_beta="",eta){
  if(is.null(groups)){
    mw=V(getMetaweb(gList))$name
    names(mw)=mw
  }else{
    mw=groups
  }

  a=disPairwise(gList, mw, type = div,eta = eta)

  b <- matrix(0,length(gList),length(gList))
  b[lower.tri(b, diag=FALSE)] <- a
  b[upper.tri(b,diag=FALSE)]<-a
  
  betadf=data.frame(b,row.names = pt_names)
  colnames(betadf)=pt_names
  write.csv2(betadf,file_beta)  
  
  return(betadf)
}

redund_fct<-function(agg_func_occur,pt_names,eta=c(0,1,2)){
  out=lapply(as.list(eta), function(e){
    hill=apply(agg_func_occur[,pt_names],2, function(x) hill_taxa(unlist(x),q = e))
    return(hill)
  })
  names(out)=as.character(eta)
  
  rf=data.frame(do.call(what = cbind,args = out))
  colnames(rf)=as.character(eta)
  return(rf)
}

#allrf=redund_fct(agg_occur,pt_names,eta=c(0,1,2))


#################### Plot functions ################################################
###PCOA
pcoa_plot=function(D,Y=NULL,corr="none",groups=NULL,pt_names=NULL,title="PCoA"){
  x=pcoa(D, correction=corr)
  data=data.frame(x$vectors[,1:2])
  colnames(data)<-c("PC1","PC2")
  data$group=groups

  p<- ggplot(data,aes(x=PC1,y=PC2,colour=group,label = row.names(data))) + 
    geom_point(size =5) +
    geom_text(col = 'black') +
    ggtitle(title)
  return(list(pcoa=x,plot=p))
}
