library(igraph)
library(data.tree)
library(miscTools)
library(parallel)

########### Knowledge graph ###################

###0. Taxonomy 
names<-read.csv("0.taxonomy.csv")
names$occurrenceId<-1:dim(names)[1]
names$key=as.character(names$key)

taxo_levels<-c("kingdom","phylum","class","order","family","genus","species")

names_tree<-names[append(taxo_levels,c("key","occurrenceId","verbatimScientificName"))]
names_tree$pathString <- paste("root2", 
                               names$kingdom, 
                               names$phylum,
                               names$class,
                               names$order,
                               names$family,
                               names$genus,
                               names$species,
                               sep = "/")

names_tree$key=as.character(names_tree$key)

###1. Trophic properties
troph<-read.csv("1.prior_interactions.csv",stringsAsFactors = FALSE)

# ###Add names to taxa
# nm=names[c("key","verbatimScientificName")]
# nm_cons=troph[c("consumer_gbif_key","consumer_name")]
# colnames(nm_cons)=c("key","verbatimScientificName")
# nm_res=troph[c("resource_gbif_key","resource_name")]
# colnames(nm_res)=c("key","verbatimScientificName")
# nm_troph=rbind(nm_cons,nm_res)
# 
# ###Check for conflictual associations
key_names=names[c("key","verbatimScientificName","rank")]

###2. µhabitat
agg<-function(x){
  if(is.factor(x)){
    return(x[1])
  }else
    return(any(x))
}
micro_hab<-read.csv("3.prior_micro_hab.csv")
micro_hab_merg<-unique(merge(key_names,micro_hab,by.x="verbatimScientificName",by.y="scientific_name","join"))
micro_habit<-aggregate(micro_hab_merg[,c(1,3:6)],by=list(key=micro_hab_merg[,2]),FUN = agg)

###3. Filling strict consumer and resource list
get_out<-function(k,type="eats"){
  res=subset(troph,consumer_gbif_key==k & interaction_type==type)$resource_gbif_key
  return(res[!is.na(res)])

  }

get_in<-function(k,type="eats"){
  res=subset(troph,resource_gbif_key==k & interaction_type==type)$consumer_gbif_key
  return(res[!is.na(res)])
}

get_microhab<-function(k){
  res=array(unlist(subset(micro_habit,key==k)[c("µhab_surf","µhab_subsurf","µhab_soil")]))
  if(length(res)==0) res=c(0,0,0)
  if(is.null(res)) res=c(0,0,0)
  return(as.integer(res))
}

##Assign to each taxa list of consumers and parasites + list of resources and hosts
names_tree$resources<-lapply(names_tree$key,get_out,"eats")
names_tree$hosts<-lapply(names_tree$key,get_out,"parasiteOf")
names_tree$consumers<-lapply(names_tree$key,get_in,"eats")
names_tree$parasites<-lapply(names_tree$key,get_in,"parasiteOf")
#names_tree$microhab<-lapply(names_tree$key,get_microhab)

###4. Build tree
taxonomy<-as.Node(names_tree)
###5. Complete µhabitat on ancestors as a logical OR of their children's µhabitat preference values
###First, we init leaves µhabitat (if not known, kept at logical(0))
###Then, we genralize the info on children to parents + info known on parents
trav=Traverse(node = taxonomy,traversal = "post-order")

UpMicroHab<-function(node){
  if(isLeaf(node)){
    ###ADDED on 17/01, on leaf inheriting global knowledge about their parent (before merge of cousins)
    nodemh=get_microhab(node$key)
    parentmh=get_microhab(node$parent$key)
    node$microhab=as.integer(nodemh|parentmh)
  }else{
    known_hab=get_microhab(node$key)
    child_microhab=lapply(node$children, function(x) res=x$microhab)
    reduced_child_mh=Reduce("|",child_microhab)
    
    node$microhab=as.integer(known_hab | reduced_child_mh)
  }  
}

#mhdf=matrix(ncol = 3, nrow = 0)
for(node in trav){
  node$microhab=UpMicroHab(node)
  # if(!is.null(node$key)){
  #   #if(node$name=="Bacteria") print("here i am")
  #   mhdf=insertRow(m=mhdf,r=nrow(mhdf)+1,v=node$microhab,rName = as.character(node$key))
  # }
}

###We update unknown microhabitats with info from closest parent
taxonomy$Do(function(node){
  node$microhab=node$parent$microhab
  #if(!is.null(node$key)) mhdf<<-insertRow(mhdf,nrow(mhdf)+1,node$microhab,rName=as.character(node$key))
  },
  filterFun = function(x) all(x$microhab==c(0,0,0)),traversal = "pre-order")

#mhdf_agg=rowsum(mhdf,row.names(mhdf))>0
print(taxonomy,"key","microhab",limit=40)

n=taxonomy$Get("name")
k=as.character(taxonomy$Get("key"))
habs=data.frame(t(taxonomy$Get("microhab")))
colnames(habs)=c("surf","subsurf","soil")
rownames(habs)<-1:nrow(habs)
habs$verbatimScientificName=n
habs$key=k

###6. Complete  trophic lists with info inherited from parents
Complete_troph<-function(node){  
  ##Traverses the taxonomic tree and on each node, appends to trophic columns the values of all their ancestors
  node$res_full=array(unlist(node$Get('resources', traversal = "ancestor")))
  node$res_full=node$res_full[!is.na(node$res_full)]
  node$cons_full=array(unlist(node$Get('consumers', traversal = "ancestor")))
  node$cons_full=node$cons_full[!is.na(node$cons_full)]
  node$host_full=array(unlist(node$Get('hosts', traversal = "ancestor")))
  node$host_full=node$host_full[!is.na(node$host_full)]
  node$paras_full=array(unlist(node$Get('parasites', traversal = "ancestor")))
  node$paras_full=node$paras_full[!is.na(node$paras_full)]
  
  return(node)
}

taxonomy$Do(Complete_troph, filterFun = isNotRoot)
print(taxonomy,"key","host_full","paras_full","res_full","cons_full",limit=20)


###7. Generate edge list
##Traversal 
Relev_node<-function(x){
  if(isNotRoot(x)){
    k=x$key
    if(is.null(k)){
      return(FALSE)
    }else{
      if(is.na(k)){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  }else{
    return(FALSE)
  }
}

Relev_node_up<-function(x){
  if(isNotLeaf(x)){
    k=x$key
    if(is.null(k)){
      return(FALSE)
    }else{
      if(is.na(k)){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  }else{
    return(FALSE)
  }
}
nodes=Traverse(taxonomy,filterFun = Relev_node)

Generate_Edges<-function(nodes){
  df <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("consumer", "resource", "type")
  colnames(df) <- x
  
  for(node in nodes){
    res=node$res_full
    cons=rep(node$key,length(res))
    type=rep("eats",length(res))
    df=rbind(df,data.frame("consumer"=cons,"resource"=res,"type"=type))
    
    cons=node$cons_full
    res=rep(node$key,length(cons))
    type=rep("eats",length(cons))
    df=rbind(df,data.frame("consumer"=cons,"resource"=res,"type"=type))
    
    host=node$host_full
    paras=rep(node$key,length(host))
    type=rep("parasiteOf",length(host))
    df=rbind(df,data.frame("consumer"=paras,"resource"=host,"type"=type)) 
    
    paras=node$paras_full
    host=rep(node$key,length(paras))
    type=rep("parasiteOf",length(paras))
    df=rbind(df,data.frame("consumer"=paras,"resource"=host,"type"=type))     
  }
  
  return(df)
}

edf=Generate_Edges(nodes)
edf=unique(edf)

write.csv2(edf,"OUT.raw_edges.csv")  ### Backup before testing parallel computations

Get_Name<-function(k,key_names){
  opt=subset(key_names,key==k)$verbatimScientificName
  # if(length(opt)>1){
  #   print(k)
  #   print("Warning possible conflicting keys")
  # }
  name=as.character(opt[1])
  return(name)
}

Get_Rank<-function(k,key_names){
  opt=subset(key_names,key==k)$rank
  # if(length(opt)>1){
  #   print(k)
  #   print("Warning possible conflicting keys")
  # }
  rank=as.character(opt[1])
  return(rank)
}

install.packages("doParallel")
install.packages("foreach")
install.packages("iterators")

library(doParallel)
numCores=detectCores()
cl <- makeCluster(numCores, type='PSOCK')
registerDoParallel(cl)

edf$consumer_name=parLapply(cl,edf$consumer,Get_Name,key_names)
edf$resource_name=parLapply(cl,edf$resource,Get_Name,key_names)

edf$consumer_rank=parLapply(cl,edf$consumer,Get_Rank,key_names)
edf$resource_rank=parLapply(cl,edf$resource,Get_Rank,key_names)

# turn parallel processing off and run sequentially again:
registerDoSEQ()

### Save edge list without habitat filtering
edf_str=apply(edf,2,as.character)
write.csv2(edf_str,"Out.interactions_inferred.csv")

# keys=as.character(rownames(mhdf_agg))
# habs=data.frame(mhdf_agg,row.names = keys)
# colnames(habs)=c("surf","subsurf","soil")
# habs$key=row.names(habs)

habs=habs[!is.na(habs$key),]
rownames(habs)=habs$key
edf$res_surf=as.numeric(habs[as.character(edf$resource),1])
edf$res_subsurf=as.numeric(habs[as.character(edf$resource),2])
edf$res_soil=as.numeric(habs[as.character(edf$resource),3])

edf$cons_surf=as.numeric(habs[as.character(edf$consumer),1])
edf$cons_subsurf=as.numeric(habs[as.character(edf$consumer),2])
edf$cons_soil=as.numeric(habs[as.character(edf$consumer),3])

#saveRDS(taxonomy,"knowledge_graph")

###Filter interactions by µhabitat
install.packages("pracma")
library(pracma)
edf$cooccur=rowSums(as.matrix(edf[,c("res_surf","res_subsurf","res_soil")])*as.matrix(edf[,c("cons_surf","cons_subsurf","cons_soil")]))/3

saveRDS(edf,"full_interactions_habit")

######Check if some habitats are uknown
nohb_res=data.frame(unique(subset(edf,is.na(res_surf))[c("resource","resource_name")]))
nohb_cons=data.frame(unique(subset(edf,is.na(cons_surf))[c("consumer","consumer_name")]))

colnames(nohb_cons)=colnames(nohb_res)=c("key","name")
nohb=rbind(nohb_cons,nohb_res)

