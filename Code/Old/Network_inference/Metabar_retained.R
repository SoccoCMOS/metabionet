
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(this.dir,"Code/Network_inference/Taxa_filtering.R",sep="/"))

error_msgs=c("Taxonomic levels are not supported","Taxonomic levels do not exist in metabar output")

retained_metabar<-function(mbraw,high="species",low="species",suffix="_name"){
  checkup=check_taxa_params(colnames_mbraw=colnames(mbraw),high=high,low=low,suffix=suffix)
  if(checkup$error==0){
    om=occurr_retained(mbraw,checkup$targetcols)
    mb=om$occur_ret
  }
  else{
    print("Wrong parameters")   
    print(error_msgs[checkup$error])
    mb=NULL
  }
  
  return(metabar=mb)
}
