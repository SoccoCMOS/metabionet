#### R6Class that represents a taxon ######
Taxon <- R6Class("Taxon", list(
  ### Attributes ###
  gbif_id = "0",
  label = "tax",
  rank = factor(c("Kingdom", "Phyllum", "Class", "Order","Family","Genus","species","undetermined")), 
  broadtype ="undetermined",
  taxref = NA,   ### Reference to node in the taxonomic tree (to get its taxonomy)
  microhab=list(surf=F,subsurf=F,soil=F),
  
  ### Object functions ###
  initialize = function(gbif_id="",label = "",rank="undetermined",broad="undetermined",microhab=list(surf=F,subsurf=F,soil=F)) {
    self$gbif_id <- gbif_id
    self$label <- label
    self$rank=rank
    self$broadtype=broad
    self$microhab=microhab
  },
  
  print = function(){
    cat("Taxon: \n")
    cat("  GBIF Key:  ", self$gbif_id, "\n", sep = "")
    cat("  Name: ", self$label, "\n", sep = "")
    cat("  Rank: ", self$rank, "\n", sep = "")
    cat(" Broad: ",self$broadtype, "\n", sep="")
    invisible(self)
  }
  ### Other functions needed ###
  ## TODO: 
  ## Function that sets the taxref using gbif_key and/or name and takes a taxonomy tree object as argument
  
))


############## Unit test ######################
taxon_object_ut<-function(){
  spec=Taxon$new("123","Species")
  spec$print()
}

taxon_object_ut()
