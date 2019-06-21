############################################################
# R script 2 of: Birds exposed to physiological stress     #
# post-breeding engage in stress-reducing social           #
# interactions in winter flocks                            #

# Running the main analysis and making the figures 3, 4, 5 #
# author: Lionel Hertzog                                   #
############################################################

# set wd to the local position of the github repo
setwd("~/Documents/PostDoc_Ghent/Feeder_stuff/network_stress/")

# load packages
library(tidyverse)
library(gridExtra)

# load data
network_dat <- read.csv("data/network_dat.csv")

# 1. Fit the centrality ~ cort_orig models
## some data points were dropped to ensure that
## model assumptions were met
m_1 <- lm(betw ~ cort_orig + sex + age + barl_orig, network_dat, subset = cort_orig < 0.1) # males have higher betw slight effect of cort
m_2 <- lm(degree ~ cort_orig + sex + age + barl_orig, network_dat[-153,]) # effect of cort
m_3 <- lm(eigen ~ cort_orig + sex + age + barl_orig, network_dat[-153,]) # trend for positive effect of cort on centrality
m_4 <- lm(strength ~ cort_orig + sex + age + barl_orig, network_dat[-153,]) # effect of cort

## save these models
m_orig <- list(betw = m_1, degree = m_2, eigen = m_3, strength = m_4)
saveRDS(m_orig, "models/m_orig.rds")

## make the figure
### the data frame to derive the prediction from
newdat <- expand.grid(cort_orig = seq(0,0.15,length=10), sex = "f", age = "1Y",barl_orig = 14.26)
### betweenness
newdat$betw <- predict(m_1, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_1)[3])) + mean(c(0,coef(m_1)[4])) # remove sex and age effect
newdat$se_betw <- predict(m_1, newdata = newdat, se.fit = TRUE)$se.fit
### degree
newdat$degree <- predict(m_2, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_2)[3])) + mean(c(0,coef(m_2)[4])) # remove sex and age effect
newdat$se_degree <- predict(m_2, newdata = newdat, se.fit = TRUE)$se.fit
### eigenvector
newdat$eigen <- predict(m_3, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_3)[3])) + mean(c(0,coef(m_3)[4])) # remove sex and age effect
newdat$se_eigen <- predict(m_3, newdata = newdat, se.fit = TRUE)$se.fit
### strength
newdat$strength <- predict(m_4, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_4)[3])) + mean(c(0,coef(m_4)[4])) # remove sex and age effect
newdat$se_strength <- predict(m_4, newdata = newdat, se.fit = TRUE)$se.fit

### the theme for all plots
theme_gg <- theme(panel.background = element_blank(), axis.line = element_line())

### the plots
gg_betw <- ggplot(newdat, aes(x=cort_orig,y=betw)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=betw - 1.96*se_betw, ymax = betw + 1.96*se_betw)) +
  geom_point(data=network_dat[-153,]) +
  labs(x="", y="Betweenness centrality") +
  theme_gg

gg_degree <- ggplot(newdat, aes(x=cort_orig,y=degree)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=degree - 1.96*se_degree, ymax = degree + 1.96*se_degree)) +
  geom_point(data=network_dat[-153,]) +
  labs(x="", y="Degree") +
  theme_gg

gg_eigen <- ggplot(newdat, aes(x=cort_orig,y=eigen)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=eigen - 1.96*se_eigen, ymax = eigen + 1.96*se_eigen)) +
  geom_point(data=network_dat[-153,]) +
  labs(x="", y="Eigenvector centrality") +
  theme_gg

gg_strength <- ggplot(newdat, aes(x=cort_orig,y=strength)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=strength - 1.96*se_strength, ymax = strength + 1.96*se_strength)) +
  geom_point(data=network_dat[-153,]) +
  labs(x="", y="Strength") +
  theme_gg

### put this all together
gg_orig <- grid.arrange(gg_eigen,gg_betw,gg_degree,gg_strength,bottom="Original feather CORT (\U003BCg)")
ggsave("figures/fig_03.png",gg_orig)

# 2. Fit the cort_induced ~ centrality models
## some data points were dropped to ensure that
## model assumptions were met
m_5 <- glm(cort_ind + 0.0001 ~ betw + sex + age + barl_ind, network_dat[-272,], family = Gamma(link="log"))
m_6 <- glm(cort_ind + 0.0001 ~ degree + sex + age + barl_ind, network_dat, subset = cort_ind < 8, family = Gamma(link="log"))
m_7 <- glm(cort_ind + 0.0001 ~ eigen + sex + age + barl_ind, network_dat[-c(25,272,279),], family = Gamma(link="log"))
m_8 <- glm(cort_ind + 0.0001 ~ strength + sex + age + barl_ind, network_dat[-272,], family = Gamma(link="log"))

## save the models
m_ind <- list(betw = m_5, degree = m_6, eigen = m_7, strength = m_8)
saveRDS(m_ind, "models/m_ind.rds")

## plot the predicted responses
### betw
newdat <- expand.grid(betw = seq(0,8.6,length=10), sex = "f", age = "AD",barl_ind = 11.55)
newdat$cort_ind <- predict(m_ind$betw, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_ind$betw)[3])) + mean(c(0,coef(m_ind$betw)[4]))
newdat$se <- predict(m_ind$betw, newdata = newdat, se.fit = TRUE)$se.fit
newdat$LCI <- exp(newdat$cort_ind - 1.96 * newdat$se)
newdat$UCI <- exp(newdat$cort_ind + 1.96 * newdat$se)
newdat$cort_ind <- exp(newdat$cort_ind)

gg_betw <- ggplot(newdat, aes(x=betw,y=cort_ind)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=LCI, ymax = UCI)) +
  geom_point(data=network_dat[-272,]) +
  labs(x="Betwenness centrality", y="Induced feather CORT (\U003BCg)") +
  ylim(c(0,2)) +
  theme_gg

### degree
newdat <- expand.grid(degree = seq(42,139,length=10), sex = "f", age = "AD",barl_ind = 11.55)
newdat$cort_ind <- predict(m_ind$degree, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_ind$degree)[3])) + mean(c(0,coef(m_ind$degree)[4]))
newdat$se <- predict(m_ind$degree, newdata = newdat, se.fit = TRUE)$se.fit
newdat$LCI <- exp(newdat$cort_ind - 1.96 * newdat$se)
newdat$UCI <- exp(newdat$cort_ind + 1.96 * newdat$se)
newdat$cort_ind <- exp(newdat$cort_ind)

gg_degree <- ggplot(newdat, aes(x=degree,y=cort_ind)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=LCI, ymax = UCI)) +
  geom_point(data=network_dat[-272,]) +
  labs(x="Degree", y="Induced feather CORT (\U003BCg)") +
  xlim(c(40,140)) +
  ylim(c(0,2)) +
  theme_gg

### eigenvector
newdat <- expand.grid(eigen = seq(0,1,length=10), sex = "f", age = "AD",barl_ind = 11.55)
newdat$cort_ind <- predict(m_ind$eigen, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_ind$eigen)[3])) + mean(c(0,coef(m_ind$eigen)[4]))
newdat$se <- predict(m_ind$eigen, newdata = newdat, se.fit = TRUE)$se.fit
newdat$LCI <- exp(newdat$cort_ind - 1.96 * newdat$se)
newdat$UCI <- exp(newdat$cort_ind + 1.96 * newdat$se)
newdat$cort_ind <- exp(newdat$cort_ind)

gg_eigen <- ggplot(newdat, aes(x=eigen,y=cort_ind)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=LCI, ymax = UCI)) +
  geom_point(data=network_dat[-c(25,272,279),]) +
  labs(x="Eigenvector centrality", y="Induced feather CORT (\U003BCg)") +
  ylim(c(0,2)) +
  theme_gg

### strength
newdat <- expand.grid(strength = seq(0,19,length=10), sex = "f", age = "AD",barl_ind = 11.55)
newdat$cort_ind <- predict(m_ind$strength, newdata = newdat, se.fit = TRUE)$fit + mean(c(0,coef(m_ind$strength)[3])) + mean(c(0,coef(m_ind$strength)[4]))
newdat$se <- predict(m_ind$strength, newdata = newdat, se.fit = TRUE)$se.fit
newdat$LCI <- exp(newdat$cort_ind - 1.96 * newdat$se)
newdat$UCI <- exp(newdat$cort_ind + 1.96 * newdat$se)
newdat$cort_ind <- exp(newdat$cort_ind)

gg_strength <- ggplot(newdat, aes(x=strength,y=cort_ind)) +
  geom_line() +
  geom_ribbon(alpha = 0.1,aes(ymin=LCI, ymax = UCI)) +
  geom_point(data=network_dat[-c(272),]) +
  labs(x="Strength", y="Induced feather CORT (\U003BCg)") +
  ylim(c(0,2)) +
  theme_gg

gg_ind <- grid.arrange(gg_eigen,gg_betw,gg_degree,gg_strength)
ggsave("figures/fig_04.png",gg_ind)

# 3. Fit the feeding frequency ~ network centrality models
m_9 <- lm(feeding_freq ~ betw + sex + age, network_dat)
m_10 <- lm(feeding_freq ~ degree + sex + age, network_dat)
m_11 <- lm(feeding_freq ~ eigen + sex + age, network_dat)
m_12 <- lm(feeding_freq ~ strength + sex + age, network_dat)

## save these models
m_feed <- list(betw = m_9, degree = m_10, eigen = m_11, strength = m_12)
saveRDS(m_feed, "models/m_feed.rds")

## plot (also with prediction interval)
### betwenness
newdat <- data.frame(betw = seq(min(network_dat$betw,na.rm=TRUE),max(network_dat$betw,na.rm=TRUE),length.out = 10), age = "1Y", sex = "f")
pred <- predict(m_9,newdat,interval = "confidence")
newdat <- cbind(newdat, pred)
names(newdat)[4] <- "feeding_freq"
pred_2 <- predict(m_9,newdat,interval="prediction")
colnames(pred_2)[2:3] <- c("P_lwr","P_upr")
newdat <- cbind(newdat, pred_2[,2:3])
cst <- mean(c(0,coef(m_9)[3])) + mean(c(0,coef(m_9)[4])) # constant to average out sex and age effects
newdat[,4:8] <- newdat[,4:8] + cst

gg_betw <- ggplot(network_dat,aes(x=betw,y=feeding_freq)) +
  geom_point() +
  geom_ribbon(data=newdat,aes(ymin=P_lwr,ymax=P_upr),color=NA,alpha=0.08) + # prediction interval
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),color=NA,alpha=0.2) + # confidence interval
  geom_line(data=newdat) +
  labs(x="Betweenness centrality", y = "Feeder visits") +
  theme_gg

### eigenvector
newdat <- data.frame(eigen = seq(min(network_dat$eigen,na.rm=TRUE),max(network_dat$eigen,na.rm=TRUE),length.out = 10), age = "1Y", sex = "f")
pred <- predict(m_11,newdat,interval = "confidence")
newdat <- cbind(newdat, pred)
names(newdat)[4] <- "feeding_freq"
pred_2 <- predict(m_11,newdat,interval="prediction")
colnames(pred_2)[2:3] <- c("P_lwr","P_upr")
newdat <- cbind(newdat, pred_2[,2:3])
cst <- mean(c(0,coef(m_11)[3])) + mean(c(0,coef(m_11)[4])) # constant to average out sex and age effects
newdat[,4:8] <- newdat[,4:8] + cst

gg_eigen <- ggplot(network_dat,aes(x=eigen,y=feeding_freq)) +
  geom_point() +
  geom_ribbon(data=newdat,aes(ymin=P_lwr,ymax=P_upr),color=NA,alpha=0.08) + # prediction interval
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),color=NA,alpha=0.2) + # confidence interval
  geom_line(data=newdat) +
  labs(x="Eigenvector centrality", y = "Feeder visits") +
  theme_gg

### degree
newdat <- data.frame(degree = seq(min(network_dat$degree,na.rm=TRUE),max(network_dat$degree,na.rm=TRUE),length.out = 10), age = "1Y", sex = "f")
pred <- predict(m_10,newdat,interval = "confidence")
newdat <- cbind(newdat, pred)
names(newdat)[4] <- "feeding_freq"
pred_2 <- predict(m_10,newdat,interval="prediction")
colnames(pred_2)[2:3] <- c("P_lwr","P_upr")
newdat <- cbind(newdat, pred_2[,2:3])
cst <- mean(c(0,coef(m_10)[3])) + mean(c(0,coef(m_10)[4])) # constant to average out sex and age effects
newdat[,4:8] <- newdat[,4:8] + cst

gg_degree <- ggplot(network_dat,aes(x=degree,y=feeding_freq)) +
  geom_point() +
  geom_ribbon(data=newdat,aes(ymin=P_lwr,ymax=P_upr),color=NA,alpha=0.08) + # prediction interval
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),color=NA,alpha=0.2) + # confidence interval
  geom_line(data=newdat) +
  labs(x="Degree", y = "Feeder visits") +
  theme_gg

### strength
newdat <- data.frame(strength = seq(min(network_dat$strength,na.rm=TRUE),max(network_dat$strength,na.rm=TRUE),length.out = 10), age = "1Y", sex = "f")
pred <- predict(m_12,newdat,interval = "confidence")
newdat <- cbind(newdat, pred)
names(newdat)[4] <- "feeding_freq"
pred_2 <- predict(m_12,newdat,interval="prediction")
colnames(pred_2)[2:3] <- c("P_lwr","P_upr")
newdat <- cbind(newdat, pred_2[,2:3])
cst <- mean(c(0,coef(m_12)[3])) + mean(c(0,coef(m_12)[4])) # constant to average out sex and age effects
newdat[,4:8] <- newdat[,4:8] + cst

gg_strength <- ggplot(network_dat,aes(x=strength,y=feeding_freq)) +
  geom_point() +
  geom_ribbon(data=newdat,aes(ymin=P_lwr,ymax=P_upr),color=NA,alpha=0.08) + # prediction interval
  geom_ribbon(data=newdat,aes(ymin=lwr,ymax=upr),color=NA,alpha=0.2) + # confidence interval
  geom_line(data=newdat) +
  labs(x="Strength", y = "Feeder visits") +
  theme_gg

### put together and save the figure
gg_feed <- grid.arrange(gg_eigen,gg_betw,gg_degree,gg_strength)
ggsave("figures/fig_05.png",gg_feed)

# 4. fit the cort_induced ~ feeding frequency model
m_13 <- glm(cort_ind + 0.0001 ~ feeding_freq + barl_ind, network_dat,family = Gamma(link="log"), subset = cort_ind < 10) # no detectable effect



