# -----------------------------------------------------------------------------
# STEP:
#    0.0   Selection of taxa at parameters level
# -----------------------------------------------------------------------------

library(pracma)

### Possible values of high and low: phylum, class, order, family, genus, species ###
tax_levels=c("phylum","class","order", "family", "genus", "species")

################################################################################################################
                    #Checkup consistency of taxonomic levels considered with provided data #
################################################################################################################

check_taxa_params<-function(colnames_mbraw,high="species",low="species",suffix="_name"){
  ### Control integrity of input parameters ###
  error=0  ##Coherent taxonomic levels => no error
  
  ## Check if taxonomic levels exist
  high_ind=which(tax_levels==high)
  low_ind=which(tax_levels==low)
  
  if(length(high_ind)==0 || length(low_ind)==0){
    error=1   ##Invalid taxonomic level, non existent
    target_col_names=NULL
  }
  else{
    ### Prevent swapping of highest and lowest taxonomic levels
    if(high_ind<low_ind){
      tmp=high_ind
      high_ind=low_ind
      low_ind=tmp
    }
    highctr=paste(tax_levels[high_ind],suffix,sep="")
    lowctr=paste(tax_levels[low_ind],suffix,sep="")
    
    if (!((highctr %in% colnames_mbraw)&(lowctr %in% colnames_mbraw))){
      error=2   ## Not present in metabar output
      target_col_names=NULL
      
    }else {
      ### Get target columns
      target_col_names=lapply(tax_levels[high_ind:low_ind],function(x) paste(x,suffix,sep=""))
    } 
  }
  return(list(error=error,targetcols=target_col_names))
}

################################################################################################################
        # Define retained taxonomic identification and ecological preference at the level: diet, µhabitat #
################################################################################################################

occurr_retained<-function(mb_raw,target_col_names){
  target_cols=mb_raw[unlist(target_col_names)]
  
  ### Construction of retained taxa ###
  retained<-vector(mode="character",length=nrow(target_cols))
  taxlevel<-vector(mode="character",length=nrow(target_cols))
  
  todrop<-list()
  
  for (i in 1:nrow(target_cols)){
    notna=which(!is.na(target_cols[i,]))
    if (length(notna)>0) {
      retained[i] <- as.character(target_cols[i,head(notna,n=1)])
      taxlevel[i] <- as.character(target_col_names[head(notna,n=1)])
    }
    else{
      retained[i] <-NA
      taxlevel[i] <-NA
      todrop <- append(todrop,i)
    }
  }
  
  ### Filtered taxa 
  #removed<-mb_raw[unlist(todrop),]
  mb<-mb_raw[-unlist(todrop),]
  ret<-retained[-unlist(todrop)]
  taxlevel<-taxlevel[-unlist(todrop)]
  
  ### Add to original matrix
  mb$retained_tax=ret
  mb$taxlevel=taxlevel
  
  return(list(occur_ret=mb))
}


##### Unit test #############
occur_test<-function(){
  out=check_taxa_params(colnames_mbraw=colnames(metabar_raw),high="genus",low="species",suffix="_name")
  if(out$error==0){
    om=occurr_retained(metabar_raw,out$targetcols)
  }
}
