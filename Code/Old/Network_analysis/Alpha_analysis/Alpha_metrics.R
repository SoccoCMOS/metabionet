# -----------------------------------------------------------------------------
# STEP:
#    Alpha analysis of ecological networks

# -----------------------------------------------------------------------------
# LIBRARIES:
library(igraph)
library(RNewsflow)
library(reshape2)
library(Rmisc)
library(dplyr)
library(stringr)

# -----------------------------------------------------------------------------
################################################################# FUNCTIONS  ###################################################################
#################################################### BEGIN(Brown/green pathways) #############################################################
###Constants
resource_light=c("R2")
resource_som=c("R1","R3")

fungal_broad=c("fungi")
bacteria_broad=c("bacterie","archae")

phyla_opport=c()
phyla_min=c()


get_phototrophs <-function(net, res_light=resource_light){
  #target=get_roots(clust,resource_light_broad,TRUE)
  return(util_haslinkto(net,resource_light,"out"))
}

get_decomposers <- function(net, res_decomp=resource_som){
  #### get list of all decomposers <- any taxa/group consuming soil organic matter
  decomp=util_haslinkto(net,resource_som,"out")
  return(decomp)
}


############################################# Broad roots ##########################################
get_taxa_fctgrp<-function(taxa_desc,broads){
  return(taxa_desc[taxa_desc$broad %in% broads,"retained_tax"])
}

get_taxa_fctgrp_test<-function(){
  taxa_desc=unique(metabar_ret[,c("retained_tax","broad")])
  resources=get_taxa_fctgrp(taxa_desc,resources_broad)
  phototrophs=get_taxa_fctgrp(taxa_desc,phototrophs_broad)
  fungis=get_taxa_fctgrp(taxa_desc,fungal_broad)
  bacterias=get_taxa_fctgrp(taxa_desc,bacteria_broad)
}

############################################# More precise roots ##########################################
util_haslinkto<-function(net,targets,dir=c("in","out")){  ###returns lists of taxa labels that are connected to targets => example here targets are all groups of resources
  vtargets=V(net)[name %in% targets]
  l=list()
  for (v in targets){
    l=append(l,neighbors(net,v,dir)$name)
  }
  
  return(unlist(l))
}

get_fungbact_decomp<-function(net,metadata){

  decomp=get_decomposers(net)
  
  #### Get fungis and fungis groups
  all_fungis=subset(metadata,broad%in%fungal_broad)[,'key']
  
  #dec_fung=intersect(all_fungis_groups,decomp)
  dec_fung=setdiff(all_fungis,decomp)
    
  #### Get bacterias and bacteria groups
  all_bacterias=subset(metadata,broad%in%bacteria_broad)[,'key']
  #all_bacterias_groups=get_roots(clust,all_bacterias,TRUE)
  
  dec_bac=intersect(all_bacterias,decomp)
  #nondec_bac=setdiff(all_bacterias_groups,decomp)
  
  #### Unclassified decomposers - cross-checking
  #unclass=setdiff(setdiff(decomp,dec_bac),dec_fung)
  
  return(list(fung=dec_fung,bact=dec_bac))
}

get_roots <- function(clust,listtaxa,format=TRUE){
  roots_unform=unique(subset(clust,clust$retained_tax %in% listtaxa)$labels)
  if(format)
    roots=unlist(lapply(roots_unform,function(x) paste("V",x,sep="")))
  else roots=roots_unform
  return(roots)
}

pathways_stat<-function(net,sources){    ### Computes pathways statistics => provide phototroph on sources for green pathways (resources for brown respectively)
  if(vcount(net)==0){
    print("Warning empty or disconnected network")
    brg_stat=list(sum_pathsize=0,mean_pathsize=0,median_pathsize=0,sd_pathsize=0,covered_taxa=NULL,nbcovered_taxa=0)
  }
  else{
    node_sources<-V(net)[name %in% sources]
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
      d=bfs(net, root, neimode = c("in"), 
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
}
#################################################### END (Brown/green pathways computing) #############################################################


#################################################### BEGIN(alpha metrics) #############################################################
topology_metrics<-function(net,metadata){
  nb_edge <- ecount(net)
  nb_node <- vcount(net)
  connectivity <- nb_edge/(nb_node * (nb_node-1))
  density <- nb_edge/nb_node
  
  ### Detailed brown green ###
  decomposers<-get_decomposers(net)
  photot <-get_phototrophs(net)
  brg <- energypathways_detail(net=net,x=decomposers,y=photot,labx="Brown",laby="Green")
  bf <- get_fungbact_decomp(net=net,metadata = metadata)
  bact<-bf$bact
  fung<-bf$fung
  bactfung <- energypathways_detail(net=net,x=bact,y=fung,labx="Bacteria",laby="Fungal")
  
  return(list(nb_edge=nb_edge,nb_node=nb_node,connectivity=connectivity,density=density,
              brown_sources=length(decomposers),light_sources=length(photot),
              bact_sources=length(bact),fung_sources=length(fung),
              brown_green_ratio=brg$RatioXY$sum,
              brownsum=brg$X$sum,greensum=brg$Y$sum,
              bacteria_fungi_ratio=bactfung$RatioXY$sum,
              bactsum=bactfung$X$sum,fungsum=bactfung$Y$sum))
}
#################################################### END (topology metrics) #############################################################

