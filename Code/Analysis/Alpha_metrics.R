# -----------------------------------------------------------------------------
# STEP:
#    Alpha analysis of ecological networks

# -----------------------------------------------------------------------------
#source("../Analysis/Util.R")

broad_light=c("light")
broad_som=c("som")

broad_fungi=c("fungi")
broad_bacteria=c("bacterie","archae")

att="broad"
scNames_opp=c()
scNames_min=c()

#################################################### BEGIN(alpha metrics) #############################################################
topoenergetic_metrics<-function(net){#,broad_light=broad_light,broad_som=broad_som,att="broad"){
  nb_edge <- ecount(net)
  nb_node <- vcount(net)
  density <- nb_edge/(nb_node * (nb_node-1))
  #density <- nb_edge/nb_node
  
  ### Detailed brown green ###
  #print("Getting phototrophs")
  phototrophs=get_energy_roots(net,attribute=att,target_lists=broad_light,filter_attribute="",filter_values=c())
  #print("Getting decomposers")
  decomposers=get_energy_roots(net,attribute=att,target_lists=broad_som,filter_attribute="",filter_values=c())
  #print("Computing brown-green ratios")
  brown_green_stats=energypathways_detail(net,y=phototrophs,x=decomposers,laby="green",labx="brown")
  
  #print("Getting bacterias")
  bacterias=get_energy_roots(net,attribute=att,target_lists=broad_som,filter_attribute=att,filter_values=broad_bacteria)
  #print("Getting fungis")
  fungis=get_energy_roots(net,attribute=att,target_lists=broad_som,filter_attribute=att,filter_values=broad_fungi)
  #print("Computing bacteria-fungi ratios")
  bact_fungi_stats=energypathways_detail(net,y=fungis,x=bacterias,laby="fungis",labx="bacterias")
  
  return(list(nb_edge=nb_edge,nb_node=nb_node,density=density,
              brown_sources=length(decomposers),light_sources=length(phototrophs),
              bact_sources=length(bacterias),fung_sources=length(fungis),
              brown_green_ratio=brown_green_stats$RatioXY$sum,
              brownsum=brown_green_stats$X$sum,greensum=brown_green_stats$Y$sum,
              bacteria_fungi_ratio=bact_fungi_stats$RatioXY$sum,
              bactsum=bact_fungi_stats$X$sum,fungsum=bact_fungi_stats$Y$sum))
}
#################################################### END (topology metrics) #############################################################

