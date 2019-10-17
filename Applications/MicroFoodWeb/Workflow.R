############  General workflow of analysis from occurrence file and metaweb #################

meta_bionet_workflow<-function(sub_edges=NA,agg_occur=NA,taxa_metadata=NA,
                               name="PROJECT",pt_names=NA,rank_filter=c("FAMILY","GENUS","SPECIES","RESOURCE"),
                               edge_ends=c("resource","consumer"),edge_attributes=c("type","cooccur"),
                               weight_att=0,type_att=0,min_prob=0.4,min_occur=1,
                               verb=1){
  
  cat("Creating node objects for the covered ",dim(taxa_metadata)[1]," taxa \n")
  ### Create node objects for each motu with metadata (rank, microhabitat) from knoweldge base ###
  taxa_nodes=apply(taxa_metadata,1, function(x){
    microhab=list(surf=x["surf"],subsurf=x["subsurf"],soil=x["soil"])
    node=Taxon$new(gbif_id=as.character(x['key']),
                   label = as.character(x['verbatimScientificName']),
                   broad=as.character(x['broad']),
                   rank=as.character(x['rank']),
                   microhab=microhab)
  })
  
  #### Select interactions with probability of co-occurrence higher than threshold ####
  if(weight_att>0){
    sel_edges=subset(sub_edges,sub_edges[edge_attributes[weight_att]]>=min_prob)
  }else{
    sel_edges=sub_edges
  }
    
  cat("Creating metaweb of trophic/parasitic interactions \n")
  metaweb=MetaWeb$new(mwebname=name, taxa=taxa_nodes, edge_list=sel_edges[,edge_ends], edge_metadata=sel_edges[,edge_attributes])
  metaweb$print()
  
  cat("Rank filter \n")
  unsupported_ranks=setdiff(unique(taxa_metadata$rank),rank_filter)
  sel_vertices=metaweb$get_taxa_rank(ranks=rank_filter)
  del_vertices=metaweb$get_taxa_rank(ranks=unsupported_ranks)
  metaweb$graph=delete.vertices(graph = metaweb$graph,del_vertices)
  
  sel_occur=agg_occur[agg_occur$key%in%sel_vertices$name,]
  cat("Projecting metaweb to local occurrence data, with minimum significant occurrence = ",min_occur,"\n")
  subnets=lapply(pt_names, function(x){
    community_compos=as.character(sel_occur[which(sel_occur[,x]>=min_occur),"key"])
    subnetwork=metaweb$project_metanetwork(sublist_taxa_keys = community_compos,net_name = x,verbose=verb) 
    return(subnetwork)
  })
  
  cat("Computing and saving network topology and energy pathways metrics \n")
  topoenerg_metrics=lapply(subnets,topoenergetic_metrics)
  topometrics_df=do.call(rbind,lapply(topoenerg_metrics, as.data.frame))
  topometrics_df$observation_id=pt_names
  
  return(list(Metaweb=metaweb,SubNetworks=subnets,TopoEnMetrics=topometrics_df))
}


############  Analysis of shared backbone accross observation points #################
shared_backbone_analysis<-function(occur_df,exp,levels,id_col,taxa_cols,cov_th=0.8){
  ##occur_dataframe: dataframe with occurrence data in the most precise scale
  ##exp: dataframe with first column as id of the most precise observations
  ##levels: columns from exp to consider in a hierarchical order (coarsest to most precise)
  ##id_col: name of column with ids of observations in exp
  ##taxa_cols: names of columns in occur_df that contain taxa information
  ##cov_th: coverage threshold in % at what coverage percentage do we consider the taxa as omnipresent 
  
  coverage=list()
  coverage$keys=occur_df$key
  for (lev in 0:3){
    cols=levels[1:(1+lev)]
    
    ### Aggregate occurences ###
    if(length(cols)==1){
      by=list(exp[,c(cols)])
    }else{
      by=as.list(exp[,c(cols)])
    }
    groups=aggregate(exp,by=by,FUN = list)[,c(id_col,cols)]
    groups$group_id=1:dim(groups)[1]
    
    group_occur=data.frame(apply(groups,1,function(x){
      #print(x)
      range=as.vector(x[[1]])
      #print(range)
      if(length(range)>1){
        rowSums(occur_df[,range])
      }else{
        occur_df[,range]
      }
    }))
    
    new_cols=lapply(groups$group_id,function(x) paste("G_",x,sep=""))
    colnames(group_occur)=new_cols
    
    group_occur[,taxa_cols]=occur_df[,taxa_cols]
    
    #cat("Common backbone of taxa shared accross observation points at ",levels[lev+1]," scale \n")
    taxa_coverage=rowSums(group_occur[,unlist(new_cols)]>0)/length(new_cols)
    backbone=group_occur[which(taxa_coverage>=cov_th),taxa_cols]
    
    coverage[levels[lev+1]]<-list(backbone)
  }
  
  return(coverage)
}

############################ Betadiversity metrics computation ###################################