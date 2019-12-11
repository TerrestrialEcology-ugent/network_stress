############################################################
# R script 3 of: Birds exposed to physiological stress     #
# post-breeding engage in stress-reducing social           #
# interactions in winter flocks                            #

# Running the permutation analysis                         #
# note that due to the large files and large computing     #
# time needed to conduct these analysis 
# author: Lionel Hertzog                                   #
############################################################

# set wd to local location of the github repo
setwd("~/Documents/PostDoc_Ghent/Feeder_stuff/network_stress/")

# load libraries
library(plyr)
library(asnipe)
library(igraph)
library(tidyverse)
library(broom)

# 1. generate the permutations
gmm_d <- readRDS("data/gmmevent_out_2016_2017.rds")

net_orig <- get_network(gmm_d$gbi)

## permutation within days within feeder
infos <- gmm_d$metadata
gbi <- gmm_d$gbi
days <- ceiling(infos$Start / 86400) # set days as integers
days[days==0] <- 1
n_perm <- 10000
# the next command take time and potentially lots of RAM, consider running it in HPC settings
net_rand <- network_permutation(gbi, association_matrix = net_orig, permutations = n_perm,
                                locations = infos$Location, within_location = TRUE, # permutate within feeder
                                days = days, within_day = TRUE) # permutate within days 


# now derive some network values
gg_li <- apply(net_rand, 1, function(net) graph_from_adjacency_matrix(net,mode="undirected",weighted = TRUE))

net_val <- lapply(gg_li, function(net) data.frame(betw = betweenness(net, directed = FALSE),
                                                  close = closeness(net, mode = "all", normalized = TRUE),
                                                  eigen = eigen_centrality(net)$vector,
                                                  deg = degree(net, mode = "all"),
                                                  strength = strength(net),
                                                  bird_tag = colnames(net_orig),
                                                  stringsAsFactors = FALSE))

# put this in a data frame
net_df <- rbind.fill(net_val)
net_df$perm <- rep(1:n_perm, each = ncol(net_orig))

# add to these the stress / feather data
cent_perm <- read.csv("data/stress_info.csv")

net_stress <- full_join(net_df, cent_perm, by = "bird_tag")
# save this file
# write.table(net_stress, "data/net_stress_perm.csv",sep = ",", row.names = FALSE)

# 

# 2. permutation on the network centrality ~ cort_orig models
# re-load the models
m_orig <- readRDS("models/m_orig.rds")

## the betwenness model
orig_b <- data.frame(term = c("(Intercept)","cort_orig","ageAD","sexm","barl_orig"), estimate = coef(m_orig$betw))

net_stress %>%
  filter(!is.na(betw)) %>%
  group_by(perm) %>%
  do(tidy(lm(betw ~ barl_orig + age + sex + cort_orig, .[-153,]))) %>%
  left_join(orig_b, by = "term") %>%
  group_by(term) %>%
  summarise(p_perm = sum(estimate.x > estimate.y) / n_perm) %>%
  mutate(p_perm = ifelse(p_perm > 0.5, 1 - p_perm, p_perm)) -> dd_betw

## the eigen model
orig_e <- data.frame(term = c("(Intercept)","cort_orig","ageAD","sexm","barl_orig"), estimate = coef(m_orig$eigen))

net_stress %>%
  filter(!is.na(betw)) %>%
  group_by(perm) %>%
  do(tidy(lm(eigen ~ barl_orig + age + sex + cort_orig, .[-153,]))) %>%
  left_join(orig_e, by = "term") %>%
  group_by(term) %>%
  summarise(p_perm = sum(estimate.x > estimate.y) / n_perm) %>%
  mutate(p_perm = ifelse(p_perm > 0.5, 1 - p_perm, p_perm)) -> dd_eigen

## the degree model
orig_d <- data.frame(term = c("(Intercept)","cort_orig","ageAD","sexm","barl_orig"), estimate = coef(m_orig$degree))

net_stress %>%
  filter(!is.na(betw)) %>%
  group_by(perm) %>%
  do(tidy(lm(deg ~ barl_orig + age + sex + cort_orig, .[-153,]))) %>%
  left_join(orig_d, by = "term") %>%
  group_by(term) %>%
  summarise(p_perm = sum(estimate.x > estimate.y) / n_perm) %>%
  mutate(p_perm = ifelse(p_perm > 0.5, 1 - p_perm, p_perm)) -> dd_degree

## the strength model
orig_s <- data.frame(term = c("(Intercept)","cort_orig","ageAD","sexm","barl_orig"), estimate = coef(m_orig$strength))

net_stress %>%
  filter(!is.na(betw)) %>%
  group_by(perm) %>%
  do(tidy(lm(strength ~ barl_orig + age + sex + cort_orig, .[-153,]))) %>%
  left_join(orig_s, by = "term") %>%
  group_by(term) %>%
  summarise(p_perm = sum(estimate.x > estimate.y) / n_perm) %>%
  mutate(p_perm = ifelse(p_perm > 0.5, 1 - p_perm, p_perm)) -> dd_strength
