#### R6Class that represents a taxonomic tree ######
TaxonomyTree <- R6Class("TaxonomyTree", list(
  ### Attributes ###
  name = "tree of life",
  tree = NA,
  path= "",
  
  ### Object functions ###
  initialize = function(dataframe, path,name = NA) {
    self$name <- name
    self$path <- path
    self$tree <- data.tree::FromDataFrameTable(dataframe,path)
  },
  
  print = function(){
    cat("Taxonomic Tree: \n")
    cat("  Name: ", self$name, "\n", sep = "")
    
    invisible(self)
  }
  ### Other functions needed ###
  
))


############## Unit test ######################