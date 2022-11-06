
# #########################################
# #########################################
# FUNCTIONS FOR BOOTSTRAP CALIBRATION
# #########################################
# #########################################


library(rpart)

# ====================
# GENERATE SOME DATA
# ====================

rdat.MARS <- function(n, p=5, model="A")
{
	X <- NULL
	for (j in 1:p) {
		x <- runif(n)
		assign(paste("x", j, sep=""), x)
		X <- cbind(X, x)
	}
	if (model=="A") mu <- 0.1*exp(4*x1) + 4/(1+exp(-20*(x2-0.5))) + 3*x3 + 2*x4 + x5
	else if (model=="B") mu <- 10*sin(pi*x1*x2) + 20*(x3-0.5)^2 + 10*x4 + x5
	else if (model=="C") mu <- 2 + 2*sign(x1 <= 0.5)*sign(x2 <=0.5)
	else stop("The arugment model= needs to be either A or B.")
	y <- mu + rnorm(n)
	dat <- data.frame(cbind(y, X))
	names(dat) <- c("y", paste("x", 1:p, sep=""))
	return(dat)
}


# ===================================================================
# FUNCTION cart() WRAPS UP STEPS FOR OBTAINING BEST TREE WITH rpart
# YET WITH FOCUS ON THE TERMINAL NODES ONLY
# ===================================================================

cart <- function(formula, data, method="anova", control=NULL,
	size.selection=c("0SE", "1SE"), plot.it=FALSE){
	if (is.null(control)) control <- rpart.control(minsplit=20, minbucket=10, 
	  maxdepth=5, cp=0, maxcompete=0,  # NUMBER OF COMPETITIVE SPLITS
        maxsurrogate=0, usesurrogate=2, surrogatestyle=0,  	# SURROGATE SPLITS FOR MISSING DATA
        xval=10)  
	tre0 <- rpart(formula=formula, data=data, method=method, control=control); 
	if (size.selection=="0SE") {
		opt <- which.min(tre0$cptable[,"xerror"])
		best.cp <- tre0$cptable[opt, "CP"]; # print(cp.best) 
		best.tree <- prune(tre0, cp = best.cp)
	} else if (size.selection=="1SE") {
		if (plot.it) plotcp(tre0, minline = TRUE) # 1SE
		cv.error <- (tre0$cptable)[,4]
		SE1 <- min(cv.error) + ((tre0$cptable)[,5])[which.min(cv.error)]      # 1SE; CAN BE EASILY MODIFIED AS aSE FOR SOME a
		position <- min((1:length(cv.error))[cv.error <= SE1]); # print(position)
		# n.size  <- (tre0$cptable)[,2] + 1  #TREE SIZE IS ONE PLUS NUMBER OF SPLITS. 
		# best.size <- n.size[position]; # best.size
		best.cp <-  sqrt(tre0$cptable[position,1]*tre0$cptable[(position-1),1]); # print(best.cp)
		# best.cp <- tre0$cptable[position,1]; print(best.cp)
		best.tree <- prune(tre0, cp=best.cp)
	}
	else stop("The values of size.selection= must be either 0SE or 1SE")
	leaf.info <- best.tree$frame[best.tree$frame$var=="<leaf>", c(2, 4:5)]
	leaf.info$sd <- sqrt(leaf.info$dev/(leaf.info$n -1))
	# THE ROW NAMES DON'T MATCH WELL WITH TERMINAL NODES 
	n.leaf <- aggregate(dat$y, by=list(best.tree$where), FUN=length); n.leaf
	leaf.info <- cbind(node=n.leaf$Group.1, leaf.info)
	# OUTPUT
	btree.size <- NROW(leaf.info)
	list(leaf=leaf.info, btree=best.tree, cp=best.cp, size=btree.size, tree0=tre0)
}


# =========================================
# SEND A TREE DOWN A DATASET AND RECOMPUTE
# =========================================

# LEAF INFO FROM fit.cart IS EXPANDED TO INCLUDE TEST SAMPLE INFO
send.down <- function(fit.cart, data, yname="y"){
	leaf <- fit.cart$leaf
	tree <- fit.cart$btree
	node <- rpart:::pred.rpart(tree, x=rpart:::rpart.matrix(data)); 
	data$node <- node
	dat.tmp <- data[order(node), c(yname, "node")]
	leaf.test <- aggregate(dat.tmp$y, by=list(dat.tmp$node), FUN=length)
	yval.test <- aggregate(dat.tmp$y, by=list(dat.tmp$node), FUN=mean)$x
	sd.test <- aggregate(dat.tmp$y, by=list(dat.tmp$node), FUN=sd)$x
	# SUMMARIZE RESULTS
	leaf.test <- cbind(leaf.test, yval.test, sd.test)
	names(leaf.test) <- c("node", "n.test", "ybar.test", "sd.test")
	leaf.info <- merge(leaf, leaf.test, by="node", all.x = FALSE)
	return(leaf.info)
}



# ==========================================================================
# BOOTSTRAP BIAS CORRECTION (BBC) FOR SD IN REGRESSION/CLASSIFICATION TREES 
# ==========================================================================

# BBC.sd <- function()

# plot.BBC <- function()


#