############# Routines to cross-manipulate the taxonomic and interaction knowledge ###############

library(data.tree)
library(tuple)
library(sets)

### Taxonomic levels supported in order of hierarchy
levels=c("root","superregnum","regnum","phylum" ,"class","order","family","genus","species")
nblevels=length(levels) ##counting root

tree_from_taxo<-function(taxo_df,levels){
  taxo_df$pathString <- paste(taxo_df[,levels],sep = "/")
  taxo_tree=as.Node(tre)
  
  return(taxo_tree)
}

get_deepest_taxa_level<-function(tree,taxa,provided_level){  ##Give data.tree object and (taxa name as a string), provided taxonomic level is a level in reversed order of deepness
  ### Create a tree from the taxonomic data.frame for more efficiency in search
  rootnode=tree$root
  fn=FindNode(rootnode,taxa)
  lev=""
  height=provided_level
  if(is.null(fn)==FALSE){
    height=fn$height
    lev=levels[nblevels-height+1]
    
    #####TODO: special cases handling
    ###Case 1: level provided is the same as found => Do nothing
    ###Case 2: level provided is deeper than what has been found, lookup the subtree
    ###Case 3: level provided is shallower than what has been found, update belief and warn
    if(provided_level>height){  #Case 3
      print(paste("Warning: provided taxonomic level is shallower than what has been found. ",
            taxa," is of level: ",lev," but provided level is ",
            levels[nblevels-provided_level+1],sep=" "))
      conflict=3
    }
    else if(provided_level<height){ #Case 2
      print(paste("Warning: provided taxonomic level is deeper than what has been found. ",
            taxa," is of level: ",lev," but provided level is ",
            levels[nblevels-provided_level+1],sep=" "))
      conflict=2
    }
    else{
      print("Matching taxonomic levels.")
      conflict=1
    }
  }
  else{
    print(paste(taxa," not found",sep=" "))
    conflict=0
    height=-1
    lev=""
  }
  
  ## Check if 
  return(list(level_name=lev,level_height=height,subtree=fn,conf=conflict))
}



unit_test<-function(){
  taxo_df=read.csv2("../../Knowledge/taxo_tree.csv")
  tree_taxo=tree_from_taxo(taxo_df,levels)
  
  int_df=read.csv2("../../Knowledge/interactions.csv")
  consumers=int_df[,c("consumer","consumer_tax_level")]
  resources=int_df[,c("resource","resource_tax_level")]
  colnames(consumers)=colnames(resources)=c("Name","Level")
  
  taxs_list=subset(unique(rbind(consumers,resources)),is.na(taxs_list$Name)==FALSE & taxs_list$Name!="")
  
  n=dim(taxs_list)[1]
  llevels=list()
  lconflicts=list()
  lsubtrees=list()
  
  for (t in 1:n){
    tl=as.character(taxs_list[t,"Level"])
    if(is.na(tl)){
      tl=""
    }
    if(tl %in% c("","NA")){  #Bu default it's a species
      h=1
    }
    else h=get_height(tl)
    nm=taxs_list[t,"Name"]
    if(h>0){ ##Not provided: do nothing
      out=get_deepest_taxa_level(tree_taxo,taxa=nm,h)
      llevels=append(llevels,out$level_name)
      lconflicts=append(lconflicts,out$conf)
      lsubstrees=append(lsubtrees,out$subtree)
    }
  }
  
  taxs_df=data.frame(taxs_list)
  taxs_df$level=llevels
  taxs_df$conflict=lconflicts
  
  OK=subset(taxs_df,taxs_df$conflict==1)
  UNFOUND=subset(taxs_df,taxs_df$conflict==0)
  DEEPER=subset(taxs_df,taxs_df$conflict==2)
  SHALLOWER=subset(taxs_df,taxs_df$conflict==3)
  
  write.csv2(as.matrix(UNFOUND),"Buffer/no_taxo_found.csv")
}


get_height<-function(lev_name){
  idx=which(levels==lev_name)[1]
  if(length(idx)==0){
    res=0
  }
  else{
    res=nblevels-idx+1
  }
  return(res)
}


get_undocumented<-function(taxa_level_list,taxodf){
  
}