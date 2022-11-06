

###########################################
# BIAS-CORRECTION OF SD IN DECISION TREES 
# TRIAL I: 
###########################################

# -----------
# MODEL A
# -----------

rm(list=ls(all=TRUE))
source("Functions-BBC.R")

set.seed(111)
nrun <- 100
B <- 200; n <- 500; Model <- "A"; 
# positive.bias <- FALSE
positive.bias <- TRUE
p <- 5; n0 <- 10000;
OUT <- NULL
TREE <- as.list(1:nrun)
for (i in 1:nrun) {
  dat <- rdat.MARS(n=n, p=p, model=Model)
  fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
  TREE[[i]] <- fit.cart$btree
  node.0 <- rpart:::pred.rpart(fit.cart$btree, x=rpart:::rpart.matrix(dat));
  test <- rdat.MARS(n=n0, p=p, model=Model)
  info.0 <- send.down(fit.cart, data=test, yname="y"); 
  sd.un <- info.0$sd
  
  # BOOTSTRAP CORRECTION
  bias <- rep(0, length(sd.un))
  for (b in 1:B){
    print(cbind(run=i, boots=b))
    id.b <- sample(1:n, size=n, replace=TRUE) 
    dat.b <- dat[id.b,]; dat.oob <- dat[-unique(id.b),]
    fit.b <- cart(y~., data=dat.b, method="anova", size.selection="1SE", plot.it=FALSE);
    info.b <- send.down(fit.b, data=dat, yname="y")  ## SHOULD USE dat.oob?
    bias.b <-  	info.b$sd.test - info.b$sd
    if (positive.bias) bias.b <- pmax(bias.b, 0)   ### NECESSARY?
    node.b <- rpart:::pred.rpart(fit.b$btree, x=rpart:::rpart.matrix(dat)); 
    tab <- table(node.0, node.b) 
    M.prop <- prop.table(tab, 1)
    bias.b <- M.prop%*%bias.b
    bias <- bias + bias.b
  }
  bias <- bias/B
  sd.co <- sd.un + bias
  out <- cbind(tree=i, info.0,bias, sd.co)
  OUT <- rbind(OUT, out)
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c("tree", "node", "n", "dev", "ybar", "sd.uncorrected", 
                   "n.test", "ybar.test", "sd.test",
                   "bias", "sd.corrected")
head(OUT)
save(OUT, TREE, file="result-ModelA.Rdat") 



# -----------
# MODEL B
# -----------

rm(list=ls(all=TRUE))
source("Functions-BBC.R")

set.seed(111)
nrun <- 100
B <- 200; n <- 500; Model <- "B"; 
# positive.bias <- FALSE
positive.bias <- TRUE
p <- 5; n0 <- 10000;
OUT <- NULL
TREE <- as.list(1:nrun)
for (i in 1:nrun) {
  dat <- rdat.MARS(n=n, p=p, model=Model)
  fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
  TREE[[i]] <- fit.cart$btree
  node.0 <- rpart:::pred.rpart(fit.cart$btree, x=rpart:::rpart.matrix(dat));
  test <- rdat.MARS(n=n0, p=p, model=Model)
  info.0 <- send.down(fit.cart, data=test, yname="y"); 
  sd.un <- info.0$sd
  
  # BOOTSTRAP CORRECTION
  bias <- rep(0, length(sd.un))
  for (b in 1:B){
    print(cbind(run=i, boots=b))
    id.b <- sample(1:n, size=n, replace=TRUE) 
    dat.b <- dat[id.b,]; dat.oob <- dat[-unique(id.b),]
    fit.b <- cart(y~., data=dat.b, method="anova", size.selection="1SE", plot.it=FALSE);
    info.b <- send.down(fit.b, data=dat, yname="y")  ## SHOULD USE dat.oob?
    bias.b <-  	info.b$sd.test - info.b$sd
    if (positive.bias) bias.b <- pmax(bias.b, 0)   ### NECESSARY?
    node.b <- rpart:::pred.rpart(fit.b$btree, x=rpart:::rpart.matrix(dat)); 
    tab <- table(node.0, node.b) 
    M.prop <- prop.table(tab, 1)
    bias.b <- M.prop%*%bias.b
    bias <- bias + bias.b
  }
  bias <- bias/B
  sd.co <- sd.un + bias
  out <- cbind(tree=i, info.0,bias, sd.co)
  OUT <- rbind(OUT, out)
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c("tree", "node", "n", "dev", "ybar", "sd.uncorrected", 
                   "n.test", "ybar.test", "sd.test",
                   "bias", "sd.corrected")
head(OUT)
save(OUT, TREE, file="result-ModelB.Rdat") 


# ----------------------
# MODEL C (TRUE TREE)
# ----------------------

rm(list=ls(all=TRUE))
source("Functions-BBC.R")

set.seed(111)
nrun <- 100
B <- 200; n <- 500; Model <- "C"; 
# positive.bias <- FALSE
positive.bias <- TRUE
p <- 5; n0 <- 10000;
OUT <- NULL
BTREE <- as.list(1:nrun)
for (i in 1:nrun) {
  dat <- rdat.MARS(n=n, p=p, model=Model)
  fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
  BTREE[[i]] <- fit.cart$btree
  node.0 <- rpart:::pred.rpart(fit.cart$btree, x=rpart:::rpart.matrix(dat));
  test <- rdat.MARS(n=n0, p=p, model=Model)
  info.0 <- send.down(fit.cart, data=test, yname="y"); 
  sd.un <- info.0$sd
  
  # BOOTSTRAP CORRECTION
  bias <- rep(0, length(sd.un))
  for (b in 1:B){
    print(cbind(run=i, boots=b))
    id.b <- sample(1:n, size=n, replace=TRUE) 
    dat.b <- dat[id.b,]; dat.oob <- dat[-unique(id.b),]
    fit.b <- cart(y~., data=dat.b, method="anova", size.selection="1SE", plot.it=FALSE);
    info.b <- send.down(fit.b, data=dat, yname="y")  ## SHOULD USE dat.oob?
    bias.b <-  	info.b$sd.test - info.b$sd
    if (positive.bias) bias.b <- pmax(bias.b, 0)   ### NECESSARY?
    node.b <- rpart:::pred.rpart(fit.b$btree, x=rpart:::rpart.matrix(dat)); 
    tab <- table(node.0, node.b) 
    M.prop <- prop.table(tab, 1)
    bias.b <- M.prop%*%bias.b
    bias <- bias + bias.b
  }
  bias <- bias/B
  sd.co <- sd.un + bias
  out <- cbind(tree=i, info.0,bias, sd.co)
  OUT <- rbind(OUT, out)
}
OUT <- as.data.frame(OUT)
colnames(OUT) <- c("tree", "node", "n", "dev", "ybar", "sd.uncorrected", 
                   "n.test", "ybar.test", "sd.test",
                   "bias", "sd.corrected")
head(OUT)
save(OUT, BTREE, file="result-ModelC.Rdat") 












#############################
# EXPLORING THE REULSTS
#############################

rm(list=ls(all=TRUE))
library(tidyverse)

load("result-ModelA.Rdat")
# load("result-ModelB.Rdat")
# load("result-ModelC.Rdat")

ls()
names(OUT); head(OUT)
#colors<-c("sd.uncorrected"="blue4","sd.corrected"="red3")
OUT %>%
  select(tree, node, sd.uncorrected, sd.test, sd.corrected) %>%
  ggplot() +
  geom_point(aes(x=sd.test, y=sd.uncorrected, colour="sd.uncorrected"), alpha=0.1) +
  #xlab("SD based on Test Samples") + ylab("SD") +
  geom_point(aes(x=sd.test, sd.corrected, colour="sd.corrected"), alpha=0.1) +
  geom_density_2d(aes(x=sd.test, y=sd.uncorrected, colour="sd.uncorrected") ) +
  geom_density_2d(aes(x=sd.test, y=sd.corrected, color="sd.corrected") ) +
  geom_abline(slope=1, intercept = 0, colour="green4", size=1)+
labs(
  x = "SD based on Test Samples",
  y = "SD",
  title = "Plot of the Uncorrected and Bias-Corrected SD",
  color ="Legend")+ theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
  
  
  
# AVERGE TREE SIZE
Bonferroni <- TRUE
# Bonferroni <- FALSE
avg.tree.size <- mean(table(OUT$tree)); 
g <- ifelse(Bonferroni, avg.tree.size, 1)

# COVERAGE WITHOUT CONSIDEIRNG MULTIPLICITY ISSUE
conf.level <- 0.95
#alpha<-c(1:200/10000)
#alpha<-0.020
alpha <- (1-conf.level)#/g    # BONFERRONI CORRECTION
OUT %>% 
  #na.exclude() %>%
  mutate(L.naive = ybar - qnorm(1-alpha/2) * sd.uncorrected/sqrt(n), 
         U.naive = ybar + qnorm(1-alpha/2) * sd.uncorrected/sqrt(n), 
         L.BBC = ybar - qnorm(1-alpha/2) * sd.corrected/sqrt(n), 
         U.BBC = ybar + qnorm(1-alpha/2) * sd.corrected/sqrt(n), 
         L.oracle = ybar - qnorm(1-alpha/2) * sd.test/sqrt(n), 
         U.oracle = ybar + qnorm(1-alpha/2) * sd.test/sqrt(n), 
         cover.naive = sign(ybar.test >= L.naive & ybar.test <= U.naive),
         cover.BBC = sign((ybar.test >= L.BBC) & (ybar.test <= U.BBC)),
         cover.oracle = sign(ybar.test >= L.oracle & ybar.test <= U.oracle)) %>%
  # select(tree, node) %>% 
  #head(n=5) %>% 
  select(cover.naive, cover.BBC, cover.oracle) %>% 
  summarise_all(mean)
#

################################################################################
######## Further Exploration ##################
###################################################
#TREE %>%
length(TREE)
load("result-ModelA.Rdat")

source("Functions-BBC.R")
Model<-"A"
n1<-1000
p <- 5

for(i in 1:length(TREE)){
  samDat<- rdat.MARS(n=n1, p=p, model=Model)
  info.samp <- send.down(TREE[[1]], data=samDat, yname="y")
  }



#conf.level <- 0.95
#alpha <- (1-conf.level)#/g    # BONFERRONI CORRECTION
alp<-c(1:200/10000)
for(alpha in 1:200){
  print(alpha)
OUT %>% 
  #na.exclude() %>%
  mutate(L.naive = ybar - qnorm(1-alpha/2) * sd.uncorrected/sqrt(n), 
         U.naive = ybar + qnorm(1-alpha/2) * sd.uncorrected/sqrt(n), 
         L.BBC = ybar - qnorm(1-alpha/2) * sd.corrected/sqrt(n), 
         U.BBC = ybar + qnorm(1-alpha/2) * sd.corrected/sqrt(n), 
         L.oracle = ybar - qnorm(1-alpha/2) * sd.test/sqrt(n), 
         U.oracle = ybar + qnorm(1-alpha/2) * sd.test/sqrt(n), 
         cover.naive = sign(ybar.test >= L.naive & ybar.test <= U.naive),
         cover.BBC = sign((ybar.test >= L.BBC) & (ybar.test <= U.BBC)),
         cover.oracle = sign(ybar.test >= L.oracle & ybar.test <= U.oracle)) #%>%
  # select(tree, node) %>% 
  #head(n=5) %>% 
  select(cover.naive, cover.BBC, cover.oracle) %>% 
  print(summarise_all(mean))
}

summarise_all(mean)


###########################################################
##  F-test#######
################
  
test<-var.test(c, b, alternative = "greater")  
test

#######################################################
ls()
hist(Alpha.Boots)
hist(Alpha.True)

dim(Alpha.Boots)
dim(Alpha.True)
  
