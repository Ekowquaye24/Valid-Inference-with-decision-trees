
###########################################
# BIAS-CORRECTION OF SD IN DECISION TREES
###########################################


# ==================================================
# INSPECT HOW SD WORKDS WITHIN EACH TERMINAL NODE 
# COMPARED TO TEST SAMPLE
# ==================================================

source("R-FunctionsBC.R")
n <- 500; n0 <- 10000; p <- 5
Model <- "A"
dat <- rdat.MARS(n=n, p=p, model=Model)
fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
test <- rdat.MARS(n=n0, p=p, model=Model)
out <- send.down(fit.cart, data=test, yname="y"); out

# COMPUTE DEVIANCE (SSE)
pred <- predict(fit.cart$btree, newdata=test, type ="vector") 
dev.test <- mean((as.numeric(test[,"y"]) - pred)^2)


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

# PLOT OF YBAR AND SD 
postscript(file="fig-inf-1.eps", horizontal=TRUE)
par(mfrow=c(1, 2), mar=rep(4,4))
plot(out$yval, out$ybar.test, type="p", col="blue", pch=19, 
	main="(a) Node Average", xlab="training", ylab="test")
abline(a=0, b=1, col="green", lwd=2)
plot(out$sd, out$sd.test, type="p", col="blue", pch=19,
	main="(b) Node SD", xlab="training", ylab="test")
abline(a=0, b=1, col="green", lwd=2)
dev.off()







# =============================================
# BOOTSTRAP BIAS CORRECTION ON NODE SD 
# ==============================================

set.seed(668)
B <- 2000
validation <- "OOB"
n <- 500; n0 <- 10000; p <- 5
Model <- "A"; nrun <- 1
positive.bias <- FALSE
# positive.bias <- TRUE
out <- NULL
for (i in 1:nrun){
	dat <- rdat.MARS(n=n, p=p, model=Model)
	fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
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
		fit.b <- cart(y~., data=dat.b, method="anova", size.selection="1SE", plot.it=FALSE);
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
	out <- cbind(sd.un, bias, sd.co, sd.gold)
}
out <- as.data.frame(out)
colnames(out) <- c("sd.un", "bias", "sd.co", "sd.gold")
out

par(mfrow=c(1, 2), mar=rep(4, 4))
plot(out$sd.gold, out$sd.un, ylab="uncorrected", xlab="gold", 
	type="p", col="blue", pch=19, main="(a)")
abline(a=0, b=1, col="green", lwd=2)
plot(out$sd.gold, out$sd.co, type="p", col="blue", pch=19,
	main="(b)", xlab="gold", ylab="corrected")
abline(a=0, b=1, col="green", lwd=2)


# PLOT THE RESULTS
postscript(file="fig-SD-A+.eps", horizontal=FALSE)
sd.gold <- out$sd.gold; sd.un <- out$sd.un; sd.co <- out$sd.co
par(mfrow=c(1, 1), mar=rep(4, 4))
ymax <- max(sd.un, sd.co); ymin <- min(sd.un, sd.co)
xmax <- max(sd.gold); xmin <- min(sd.gold)
plot(c(xmin, xmax), c(ymin, ymax), type="n", xlab="Test Sample SD", ylab="SD")
segments(x0=sd.gold, y0=sd.co, x1=sd.gold, y1=sd.un, lwd=.5, col="gray85")
points(sd.gold, sd.un, col="gray75", pch=16, cex=.8)
points(sd.gold, sd.co, col="gray25", pch=15, cex=.8)
legend(x=xmin, y=ymax, pch=c(16, 15), col=c("gray75", "gray25"), 
	legend=c("Uncorrected", "Bias-Corrected"), cex=1.2)
abline(a=0, b=1, col="green4", lwd=2)
dev.off()
































# ==================================================================================
# NOT WORKING --- BOOTSTRAP BIAS CORRECTION ON OVERALL VARIANCE (HOMOSCEDASTICITY)
# ==================================================================================

B <- 1000
validation <- "OOB"
n <- 500; n0 <- 10000; p <- 5
Model <- "A"; nrun <- 50
out <- NULL
for (i in 1:nrun){
	dat <- rdat.MARS(n=n, p=p, model=Model)
	fit.cart <- cart(y~., data=dat, method="anova", size.selection="1SE", plot.it=FALSE);
	leaf.info <- fit.cart$leaf
	V.un <- sum(leaf.info$dev)/(sum(leaf.info$n) - NROW(leaf.info))

	# BOOTSTRAP CORRECTION
	bias <- 0
	for (b in 1:B){
		id.b <- sample(1:n, size=n, replace=TRUE) 
		dat.b <- dat[id.b,]; dat.oob <- dat[-unique(id.b),]
		fit.b <- cart(y~., data=dat.b, method="anova", size.selection="1SE");
		info.b <- fit.b$leaf	
		V.b <- mean(info.b$dev)
		if (validation == "OOB") {
			pred.b <- predict(fit.b$btree, newdata=dat.oob, type ="vector") 
			V.0 <- mean((dat.oob$y - pred.b)^2)
		} else {
			pred.b <- predict(fit.b$btree, newdata=dat, type ="vector") 
			V.0 <- mean((dat$y - pred.b)^2)
		}
		bias <- bias + (V.0-V.b)
	}
	bias <- bias/B
	V.co <- V.un - bias

	# TEST SAMPLE BASED - GOLD
	test <- rdat.MARS(n=n0, p=p, model=Model)
	pred <- predict(fit.cart$btree, newdata=test, type ="vector") 
	V.gold <- mean((test$y - pred)^2)
	out <- rbind(out, c(i, V.un, bias, V.co, V.gold))
}
out <- as.data.frame(out)
colnames(out) <- c("run", "V.un", "bias", "V.co", "V.gold")
out



#