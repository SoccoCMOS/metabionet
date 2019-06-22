library(dplyr)
library(rlist)

int_list<-read.csv2("Knowledge/bise_int_list.csv")
taxa_list <- read.csv2("Knowledge/taxa_list.csv")

get_id_at_level<-function(t,l,k){  ##Do this for each parent of the taxa (down to its level and up to the highest relevant levels considered)
  id_lev=""
  if(!is.null(k))
    kRow<-subset(taxa_list, key==k)
  
  if(!is.null(t))
    tRow<-subset(taxa_list, verbatimScientificName==t)
  
  if(kRow$occurrenceId[1]!=tRow$occurrenceId[1])   ##TODO: check for duplicate keys and synonyms
    { print("name_key_conflict")
  }else{
    name_lev=as.character(kRow[1,l])
    key_lev=subset(taxa_list,verbatimScientificName==name_lev)$key
    if(!is.na(key_lev))
      key_lev=as.character(key_lev)
  }
  
  return(list(name=name_lev,key=key_lev))
}

get_resources_of_ids<-function(k,n,l){ #k, n and l are lists of the same size
  ###Get resources list by its key
  keyRow<-subset(int_list,consumer_gbif_key %in% k)
  
  ###Get resources list by its name
  target=data.frame(cbind(consumer_rank=toupper(l),consumer_name=n))
  nameRow<-merge(x=int_list,y=target)
  
  ###Concat both rows
  allResources<-distinct(rbind(keyRow,nameRow))
  return(allResources[c("resource_gbif_key","resource_rank","interaction_type")])
}

get_consumers_of_ids<-function(k,n,l){ #k, n and l are lists of the same size
  ###Get resources list by its key
  keyRow<-subset(int_list,resource_gbif_key %in% k)
  
  ###Get consumers list by its name
  target=data.frame(cbind(resource_rank=toupper(l),resource_name=n))
  nameRow<-merge(x=int_list,y=target)
  
  ###Concat both rows
  allConsumers<-distinct(rbind(keyRow,nameRow))
  return(allConsumers[c("consumer_gbif_key","consumer_rank","interaction_type")])
}

get_all<-function(t,k,relev_lev){  ##Returns for a taxa all names and keys of parents at relevant levels and their res/cons
  full_info<-list()
  for (l in relev_lev){
    p=get_id_at_level(t,tolower(l),k)
    res=get_resources_of_id(p$key,p$name,l)
    cons=get_consumers_of_ids(p$key,p$name,l)
    full_info[[l]]$parent=p
    full_info[[l]]$resources=res
    full_info[[l]]$consumers=cons
  }
  return(parent_l)
}                                                                                                                                                                                     