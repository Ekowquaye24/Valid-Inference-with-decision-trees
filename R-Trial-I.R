
# ====================================================
# INSPECT HOW MEAN & SD WORK WITHIN EACH TERMINAL NODE 
# COMPARED TO TEST SAMPLE
# ==================================================

rm(list=ls(all=TRUE))
source("Functions-BBC.R")
n <- 500; n0 <- 10000; p <- 5
Model <- "A"
dat <- rdat.MARS(n=n, p=p, model=Model)
fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
test <- rdat.MARS(n=n0, p=p, model=Model)
out <- send.down(fit.cart, data=test, yname="y"); out

# PLOT OF YBAR AND SD 
postscript(file="fig-inf.eps", horizontal=TRUE)
par(mfrow=c(1, 2), mar=rep(4,4))
plot(out$ybar.test, out$yval, type="p", col="blue", pch=19, 
	main="(a) Node Average", ylab="training", xlab="test")
abline(a=0, b=1, col="green", lwd=2)
plot(out$sd.test, out$sd, type="p", col="blue", pch=19,
	main="(b) Node SD", ylab="training", xlab="test")
abline(a=0, b=1, col="green", lwd=2)
dev.off()



###########################################
# BIAS-CORRECTION OF SD IN DECISION TREES 
# TRIAL I: ONE SIMULATION RUN
###########################################

rm(list=ls(all=TRUE))
source("Functions-BBC.R")

set.seed(543)
B <- 500; n <- 500; Model <- "A"; 
# positive.bias <- FALSE
positive.bias <- TRUE

p <- 5; n0 <- 10000; 
dat <- rdat.MARS(n=n, p=p, model=Model)
fit.cart <- cart(y~., data=dat, method="anova", size.selection="0SE", plot.it=FALSE);
node.0 <- rpart:::pred.rpart(fit.cart$btree, x=rpart:::rpart.matrix(dat));
test <- rdat.MARS(n=n0, p=p, model=Model)
info.0 <- send.down(fit.cart, data=test, yname="y"); 
sd.un <- info.0$sd
sd.gold <- info.0$sd.test

# BOOTSTRAP CORRECTION
bias <- rep(0, length(sd.un))
for (b in 1:B){
	print(b)
	id.b <- sample(1:n, size=n, replace=TRUE) 
	dat.b <- dat[id.b,]; dat.oob <- dat[-unique(id.b),]
	fit.b <- cart(y~., data=dat.b, method="anova", size.selection="0SE", plot.it=FALSE);
	info.b <- send.down(fit.b, data=dat, yname="y")
	bias.b <-  	info.b$sd.test - info.b$sd
	if (positive.bias) bias.b <- pmax(bias.b, 0)
	node.b <- rpart:::pred.rpart(fit.b$btree, x=rpart:::rpart.matrix(dat)); 
	tab <- table(node.0, node.b) 
	M.prop <- prop.table(tab, 1)
	bias.b <- M.prop%*%bias.b
	bias <- bias + bias.b
}
bias <- bias/B
sd.co <- sd.un + bias
out <- cbind(bias, sd.co, sd.gold)
out <- as.data.frame(cbind(fit.cart$leaf, out)); 
colnames(out) <- c("node", "n", "dev", "mean.y", "sd.un", "bias", "sd.co", "sd.gold")
out 


# PLOT THE RESULTS
# postscript(file="fig-sd-A+.eps", horizontal=FALSE)
windows()
sd.gold <- out$sd.gold; sd.un <- out$sd.un; sd.co <- out$sd.co; n.t <- out$n
scale.adj <- 150
par(mfrow=c(1, 1), mar=rep(4, 4))
ymax <- max(sd.un, sd.co); ymin <- min(sd.un, sd.co)
xmax <- max(sd.gold); xmin <- min(sd.gold)
plot(c(xmin, xmax), c(ymin, ymax), type="n", xlab="Test Sample SD", ylab="SD")
segments(x0=sd.gold, y0=sd.co, x1=sd.gold, y1=sd.un, lwd=.5, col="gray85")
points(sd.gold, sd.un, col="gray75", pch=16, cex=.8)
# text(sd.gold+ (mean(sd.gold)/scale.adj), sd.un, col="gray75", labels =n.t, cex=.5)
points(sd.gold, sd.co, col="gray25", pch=15, cex=.8)
text(sd.gold+(mean(sd.gold)/scale.adj), sd.co, col="gray25", labels=n.t, cex=.5)
legend(x=xmin, y=ymax, pch=c(16, 15), col=c("gray75", "gray25"), 
	legend=c("Uncorrected", "Bias-Corrected"), cex=1.2)
abline(a=0, b=1, col="green4", lwd=2)
# dev.off()













#