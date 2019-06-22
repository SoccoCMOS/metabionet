#### R6Class that represents a metaweb of interactions ######

MetaWeb <- R6Class("MetaWeb", list(
  ### Attributes ###
  webname = "Metaweb0",
  nodes = NA,  ### List of objects of type Taxon representing taxa identified
  graph = NA, ### Actual interaction network: igraph object
  
  ### Object functions ###
  initialize = function(mwebname, taxa, edge_list, edge_metadata=NA) {
    self$webname <- mwebname
    self$nodes <- taxa
    ### Create empty graph ###
    self$graph = graph.empty()%>%
      ### Add vertices ###
      add_vertices(nv=length(taxa))%>%
      ### Add vertex attributes ###
      set_vertex_attr("name",value=unlist(lapply(taxa,FUN=function(x) x$gbif_id)))%>%
      set_vertex_attr("scName",value=unlist(lapply(taxa,FUN=function(x) x$label)))%>%
      set_vertex_attr("rank",value=unlist(lapply(taxa,FUN=function(x) as.character(x$rank))))%>%
      set_vertex_attr("broad",value=unlist(lapply(taxa,FUN=function(x) x$broadtype)))%>%
      
      ### Add edges and their attributes ###
      add_edges(edges=as.vector(t(as.matrix(edge_list[,c("resource","consumer")]))),
                attr = list(type=as.character(edge_metadata[,"type"]),weight=as.double(edge_metadata[,"cooccur"]))
      )
  },
  
  print = function(){
    cat("Metaweb: \n")
    cat("  Name: ",self$webname, "\n", sep = "")
    cat(" Network statistics: \n")
    cat(" Number of taxa: ",igraph::vcount(self$graph),"\n", sep = "")
    cat(" Number of interactions: ",igraph::ecount(self$graph),"\n", sep = "")
    cat(" Connectivity: ",edge_density(self$graph))
    invisible(self)
  },
  
  ### Other functions needed ###
  print_network = function(){  ##TODO: add graphical attributes
    plot.igraph(graph,label="scName")
  },
  
  save_network = function(file_name,fmt= c("edgelist", "pajek", "ncol", "lgl",
                                            "graphml", "dimacs", "gml", "dot", "leda")){
    write_graph(graph=self$graph, file=file_name, format=fmt)
  },
  
  ### Add attributes to vertices ###
  add_vertex_atts=function(att_name,idx,atts){
    self$graph=set_vertex_attr(graph=self$graph,name=att_name,index=idx,value=atts)
  },
  
  ### Network projection functions => returns subnetwork of type igraph ###
  project_metanetwork=function(sublist_taxa_keys,net_name,verbose=0){  ### Project to subset of taxa provided in a list of keys in string format
    subnetwork=induced_subgraph(self$graph,vids=sublist_taxa_keys,impl="create_from_scratch")
    if(verbose==1){
      cat("Projected trophic web: \n")
      cat("  Name: ",net_name, "\n", sep = "")
      cat(" Network statistics: \n")
      cat(" Number of taxa: ",vcount(subnetwork),"\n", sep = "")
      cat(" Number of interactions: ",ecount(subnetwork),"\n", sep = "")
      cat(" Connectivity: ",edge_density(subnetwork),"\n")
    }

    return(subnetwork)
  },
  
  get_taxa_broad=function(broad_list){ ### Returns list of taxa of a certain broad type
    vbroad=V(self$graph)[broad%in%broad_list]
    return(vbroad)
  },
  
  get_taxa_rank=function(ranks){ ### Filter taxa of certain rank (provide ranks to keep)
    vrank=V(self$graph)[rank%in%ranks]
    return(vrank)
  },
  
  get_interactions_type=function(typ=c("eats","parasiteOf")){ ### Filter to a certain type of interactions
    etype=E(self$graph)[type%in%typ]
    return(etype)
  },
  
  get_interactions_strength=function(th=0){ ### Filter to a certain type of interactions
    estrength=E(self$graph)[weight>th]
    return(estrength)
  },
  
  subgraph_edges=function(eids,delete_taxa=F){
    return(subgraph.edges(self$graph, eids, delete.vertices = delete_taxa))
  },
  
  get_adjacency_matrix=function(edge_att="weight",edges=T){
    return(as_adjacency_matrix(self$graph,attr = edge_att,edges = edges,sparse = F))
  }
))
