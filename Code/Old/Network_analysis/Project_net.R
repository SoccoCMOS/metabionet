# -----------------------------------------------------------------------------
# PROJECT:
#    BISE : Linking soil foodweb and soil functions dynamics under
#    agricultural practices (des)intensification and different input levels
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# STEP:
#    Projecting metanetworks
# -----------------------------------------------------------------------------


library(RNewsflow)
library(reshape2)
library(Rmisc)
library(dplyr)
library(stringr)

###Metanetwork = adjacency matrix with group names on rows and columns
###Occur= matrix with 3 columns: retained_taxa, cluster where it is affected, m columns for the m sites with 0/1 values for respectively read/not read

project_net<-function(net,mb,idx_tax,idxs_plots,levels,clust,adj_mode=FALSE){   #net=metanerwork, mb=occurrence matrix with columns retained_taxa and point replicates, clust=labels of clustering
                                      ############### Get occurrence data at groups or mesoscopic level ##############################
  
  mb=metabar_ret[,c(idx_tax,idxs_plots)]
  ## Uniformiser les noms des lignes et colonnes (formatage)
  #row.names(net)=as.character(1:length(net)) 
  colnames(net)=row.names(net)
  ## Création de la edge list A partir de la matrice d'adjacence
  net_edges=melt(as.matrix(net))
  
  ### Projection du m?tar?seau sur la parcelle ###
  ## Filtrer les taxons qui interagissent uniquement + leurs relevés de read metabar 
  mb_short <- mb[mb$retained_tax %in% clust$retained_tax,]
  
  ## Ajouter les informations de clustering: classe attribuée
  ## => Dans mb_long: pour chaque taxon (au niveau i), on aura les reads dans les 256 points d'échantillonnages (16*4*4) + différentes classifications
  mb_long <- merge(clust, mb_short, "retained_tax")
  
  ## Aggregate number of reads by summing the value of all rows with same class/codeFW/retained_taxa
  ## mb_class = aggregated metabar reads per class
  mb_class <- aggregate(mb_long[, -c(1:2)], list(class = mb_long[,"labels"]), sum)
  
  
  ## transpose to have metabar sampling points on rows and classes on cols
  mb_class_t0 <- as.data.frame(t(mb_class[,-1]))

  ## Rename columns
  colnames(mb_class_t0) <- colnames(net)  
  
  ### mb_class_t1 [i,j] = Number of reads of MOTU classified in class j seen in point i => group per point (replicates up to point)
  ## Pour chaque plot (parcelle + traitement)

  mb_class_t1 <- aggregate(mb_class_t0, by = levels, FUN = sum)
  # 
  # #A l'echelle de la parcelle --> grouper les 4 ?chantillons ? la parcelle (points up to plots)
  mb_list_t <- split(mb_class_t1[,-1], with(mb_class_t1,level))
  mb_list <- lapply(mb_list_t, FUN = t)
  
  adj <- lapply(mb_list,mynetwork,net_edges)
  
  return(list_adj=adj)
}

mynetwork <- function(mat,meta_inter){
  step1 <- ifelse(mat>=1,1,0) ## Binarize, if at least one read of a taxa in group on row, say true else false ???? TODO: >1 or 100 or 0 ?
  step2 <- as.data.frame(step1[rowSums(step1)>0,]) ## Number of points where the considered class is read
  interactions <- meta_inter[meta_inter[,1] %in% rownames(step2) & meta_inter[,2] %in% rownames(step2),]
  filt_inter <- subset(interactions,interactions$value>0)
  if(length(filt_inter)==0){
    print("stop")
  }
  g <- graph_from_data_frame(filt_inter, directed = TRUE)
  g <- simplify(g, remove.multiple = F, remove.loops = T)
  
  adj_mat=as_adj(g)
  
  return(list(adj=adj_mat,graph=g))
}


project_net_test<-function(){
  partition$retained_tax=row.names(partition)
  project_net(net=meta_clust_adj_filt,mb=metabar_ret,idx_tax=269,idxs_plots=13:268,nbrep=4,clust=partition)
}