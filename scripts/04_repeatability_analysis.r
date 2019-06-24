############################################################
# R script 4 of: Birds exposed to physiological stress     #
# post-breeding engage in stress-reducing social           #
# interactions in winter flocks                            #

# Running the repeatability analysis                       #

# author: Lionel Hertzog                                   #
############################################################

# set wd to local position of github folder
setwd("~/Documents/PostDoc_Ghent/Feeder_stuff/network_stress/")

# load libraries
library(plyr)
library(tidyverse)
library(igraph)
library(rptR)

# load the daily networks
lif <- list.files(path = "data/daily_networks/", pattern = "out")
li_net <- list()
for(i in seq_along(lif)){
  tmp <- read.table(file = paste0("data/daily_networks/",lif[i]),head=TRUE,sep=" ",stringsAsFactors = FALSE)
  li_net[[i]] <- graph_from_adjacency_matrix(as.matrix(tmp),mode="undirected",weighted = TRUE)
  names(li_net)[i] <- substr(lif[i],8,15)
  rm(tmp)
}

#compute centrality metrics for each bird, for each day
centr <- ldply(li_net,function(g) data.frame(PIT_TAG=attr(strength(g),"names"),betw = betweenness(g, directed = FALSE),
                                             eigen = eigen_centrality(g)$vector,
                                             deg = degree(g, mode = "all"),
                                             strength = strength(g)))

#load bird informations
bird <- read.table("data/bird_infos.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
bird$Year <- ifelse(bird$Period == "Autumn2015","2016",ifelse(bird$Period == "Autumn2016","2017",""))

# some data formatting
centr %>%
  rename(Date = .id) %>%
  filter(!is.na(betw) | !is.na(strength) | !is.na(deg)) %>%
  mutate(PIT_TAG = substr(PIT_TAG,2,11)) %>%
  mutate(Year = ifelse(Date < "20160430","2016","2017")) %>%
  left_join(bird[,c(2,3,9:18)],by=c("PIT_TAG","Year")) %>%
  filter(!is.na(MASS) & !is.na(SEX) & !is.na(AGE) & !is.na(SPEC)) -> cent_dd

# add network size for every day
size_l <- ldply(li_net, function(g) data.frame(Size = vcount(g)))
size_l %>%
  rename(Date = .id) %>%
  right_join(cent_dd, by = "Date") %>%
  filter(Date > 20161201) -> cent_dd2

# compute repeatability
m_rpt_b <- rpt(betw ~ Size + (1|PIT_TAG),grname = "PIT_TAG",data = cent_dd2, datatype = "Gaussian", npermut = 1000)
m_rpt_e <- rpt(eigen ~ Size + (1|PIT_TAG),grname = "PIT_TAG",data = cent_dd2, datatype = "Gaussian", npermut = 1000)
m_rpt_d <- rpt(deg ~ Size + (1|PIT_TAG),grname = "PIT_TAG",data = cent_dd2, datatype = "Gaussian", npermut = 1000)
m_rpt_s <- rpt(strength ~ Size + (1|PIT_TAG),grname = "PIT_TAG",data = cent_dd2, datatype = "Gaussian", npermut = 1000)

# the output plot
png("figures/repeatability_permutations.png",width=1200,height=1200)
par(mfrow=c(2,2),mar=c(5.1,5.1,4.1,2.1))
plot(m_rpt_b,type="permut",main="Permutation repeatabilities for\nbetwenness, mean repeatability: 0.16 [0.13 - 0.20]",cex.axis=1.7,cex.lab=1.7,cex.main=1.8)
plot(m_rpt_e,type="permut",main="Permutation repeatabilities for\neigenvector, mean repeatability: 0.21 [0.17 - 0.25]",cex.axis=1.7,cex.lab=1.7,cex.main=1.8)
plot(m_rpt_d,type="permut",main="Permutation repeatabilities for\ndegree, mean repeatability: 0.35 [0.30 - 0.40]",cex.axis=1.7,cex.lab=1.7,cex.main=1.8)
plot(m_rpt_s,type="permut",main="Permutation repeatabilities for\nstrength, mean repeatability: 0.39 [0.34 - 0.45]",cex.axis=1.7,cex.lab=1.7,cex.main=1.8)
dev.off()

