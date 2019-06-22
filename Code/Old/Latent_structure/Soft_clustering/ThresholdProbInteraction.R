# -----------------------------------------------------------------------------
# STEP:
#    1.5    Choice of probabiliy threshold in SBM adjacency matrices
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LIBRARIES:
library(igraph)
# -----------------------------------------------------------------------------

## Criterion 1: Keep th_n ratio of nodes (trophic groups)
## Criterion 2: Keep graph connexity
## Criterion 3: No intraguild interactions

get_threshold<-function(SBM_adj,th_n=1,max_nbintra=0){  ##sqmat=square matrix of adjacency, th_n ratio of nodes to keep, max_nbinta=max number of intraguild allowed
  p_range <- seq(0,0.995,by=0.001)
  # Following metrics are computed to keep track of changes in graph entropy as we filter edges based on weight threshold p
  densities<-vector(mode="integer",length=length(p_range)) # Graph density = nb_edges/nb_nodes
  n=length(SBM_adj)
  
  optim_reached=FALSE
  node_preserv=TRUE #Criterion 1
  intrag=TRUE # criterion 3
  connex=TRUE # criterion 2
  cpt=1
  while (!optim_reached){
    cpt=cpt+1
    p=p_range[cpt]
    ## Binarize
    SBM_adj_filt <- as.data.frame(ifelse(SBM_adj > p , 1, 0))
    ## Filter
    condition <- as.data.frame(rowSums(SBM_adj_filt))
    condition$colSums <- colSums(SBM_adj_filt)
    condition$cond <- rowSums(condition)
    SBM_adj_filt <- SBM_adj_filt[condition$cond  > 0, condition$cond > 0]
    
    # Creation du graphe
    g<-graph.adjacency(as.matrix(SBM_adj_filt))
    
    # Re-compute metrics
    densities[cpt-1]=graph.density(g)
    
    # Node preservation criterion (1)
    node_preserv=(vcount(g)/n)>=th_n 
    
    # Check for intraguild criterion (3)
    intra=lapply(1:length(SBM_adj_filt),FUN=function(x) SBM_adj[x,x]>p)
    n_intrag=length(which(intra==TRUE))
    intrag=n_intrag>max_nbintra
    
    # Check for connexity criterion (2)
    clu<-components(g)
    connex=clu$no==1
    
    if(!(node_preserv & intrag & connex)){
      optim_reached=TRUE
      p=p_range[cpt-1]
    }
  }
  
  return(p)
}


get_thresholded_mat_bin<-function(mat,p,bin=TRUE){
  if(bin==FALSE){
  mat_filt<- mat
  mat_filt[mat_filt<=p] <- 0
  }
  else{
    mat_filt <- as.data.frame(ifelse(mat > p , 1, 0))
  }
  
  condition <- as.data.frame(rowSums(mat_filt))
  condition$colSums <- colSums(mat_filt)
  condition$cond <- rowSums(condition)
  
  mat_out <- mat_filt[condition$cond  > 0, condition$cond > 0]
  
  return(mat_out)
}


#########################################################

threshold_test<-function(){
  p_optim=get_threshold(meta_clust_adj,th_n=0.95,max_nbintra=0)
  mat=get_thresholded_mat_bin(meta_clust_adj,p_optim,TRUE)
}

