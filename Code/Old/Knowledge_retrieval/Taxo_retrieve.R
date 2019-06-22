############## Retrieve data from GBIF ###################
library(rgbif)

taxo_ranks=taxrank()

get_key_rank<-function(names){  ###taskes as input a list of taxa names and returns their GBIF key, their canonical name and their rank
  res=data.frame()
  unr=list()
  for(name in names){
    out = name_suggest(q=name)
    if(dim(out)[1]==0){  ###No result
      unr=append(unr,name)
    }
    else{
      l=data.frame(q=name,canonicalName=out$canonicalName[1],key=out$key[1],rank=out$rank[1])
      res=rbind(res,l)
    }
  }
  return(list(results=res,unresolved=unr))
}

get_taxo_tree<-function(key,rnk){   ##TODO: cross-check this one for particular cases
  hier=occ_search(taxonKey=key, return="hier")[[1]]
  rank_nb=which(taxo_ranks==tolower(rnk))
  
  return(hier[1:rank_nb,])
}

unit_test<-function(){
  names=c("lumbricidae","pterostichus","xhekjgheg","bird")
  out1=get_key_rank(names)$results
  
  for(i in 1:dim(out1)[1]){
    hieri=get_taxo_tree(out1$key[i],out1$rank[i])
  }
}