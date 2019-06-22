############  Installing dependencies #################
list.of.packages <- c("ggplot2", "hillR","entropy","econetwork","Rcpp","R6","data.tree","igraph","blockmodels","RNewsflow","reshape2","Rmisc","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

############  Load libraries #################
library(R6)
library(data.tree)
library(igraph)
library(RNewsflow)
library(reshape2)
library(Rmisc)
library(dplyr)
library(stringr)
library(blockmodels)
library(econetwork)
library(entropy)
library(hillR)
library(ape)


##### Source all needed files #####
source("Code/Data Structures/Taxon.R")
source("Code/Data Structures/Taxonomy.R")
source("Code/Data Structures/MetaWeb.R")
source("Code/Analysis/Util.R")
source("Code/Analysis/Alpha_metrics.R")
source("Applications/BISE/Workflow.R")