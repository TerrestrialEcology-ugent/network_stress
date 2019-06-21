##########################################################
# R script 1 of: Birds exposed to physiological stress   #
# post-breeding engage in stress-reducing social         #
# interactions in winter flocks                          #

# Generating social network and deriving network metric  #

# author: Lionel Hertzog                                 #
##########################################################

# set wd to the local position of the github repo
setwd("~/Documents/PostDoc_Ghent/Feeder_stuff/network_stress/")

# load packages
library(asnipe)
library(igraph)
library(reshape2)

# 1. Generate the network from the raw feeder data
## note: this analysis takes quite some RAM and some time
## it is included here only for reproducability, you can
## directly get the network from the data/network.csv file

# # load data
feeder_dat <- read.csv("data/feeder_data_formatted.csv")

# get the infos for computing the social network
## the time of the recording
day_hour <- strptime(feeder_dat$day_hour, format = "%Y-%m-%d %H:%M:%S")
day_hour <- as.numeric(day_hour) # turn into integer
day_hour <- day_hour - min(day_hour) # for faster computation

## the transponder id of the bird
bird_id <- feeder_dat$Transponder.code

## the feeder id
feeder_id <- feeder_dat$Feeder_ID

# run gmmevents, this line of code takes a lot of RAM! Beware !!!
# gmm <- gmmevents(day_hour, bird_id, feeder_id, verbose = FALSE)

# get network
# network <- get_network(gmm$gbi)

# # save it in a file
# # write.table(network, "data/network.csv",row.names = FALSE, sep = ",")

# 2. Derive the network indices

## load the network
network <- read.csv("data/network.csv")

## turn it into an igraph object
gg_winter <- graph_from_adjacency_matrix(as.matrix(network), mode = "undirected", weighted = TRUE)

## put all into one data frame
network_metric <- data.frame(bird_tag = gsub("X","",names(V(gg_winter))),
                             betw = log(betweenness(gg_winter, directed = FALSE) + 1),
                             eigen = eigen_centrality(gg_winter)$vector,
                             degree = degree(gg_winter, mode = "all"),
                             strength = strength(gg_winter),
                             stringsAsFactors = FALSE)

# 3. Merge with the stress and other covariate

## load the stress info
stress <- read.csv("data/stress_info.csv")

## merge
network_dat <- merge(network_metric, stress, by = "bird_tag")

## save it in a file
write.table(network_dat,"data/network_dat.csv", row.names = FALSE, sep = ",")


