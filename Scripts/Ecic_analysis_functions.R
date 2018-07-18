# Functions for Erodium CA/Chile paper
# Andrew Latimer April 2017

# function to plot results of anova as a stacked bar chart
aov.stackbar <- function(a, a.names) { # a is a list of anova models
  n = length(a)
  nfactors = length(anova(a[[1]])[,1])-1
  dtemp = matrix(0, nfactors, n)
  for (i in 1:n) {
    z = anova(a[[i]])
    SST = sum(z[,2])
    SS = z[,2]/SST
    dtemp[,i] = SS[1:nfactors]
  }
  par(mar=rep(5, 4), pty="s")
  plot.new()
  pal <- wes_palette("FantasticFox1")
  barplot(dtemp, ylim=c(0,1), col=pal,  names.arg=a.names, cex.names=1.1, ylab="Proportion variance explained", cex.axis=1.2, cex.lab=1.4, legend.text=c("Population", "Treatment", "Pop x Trmt"), args.legend=list(x="topright"))
}


# plot random effects of a model with a single set of random intercepts
ranefplot <- function(m, axis.label) {
  randoms<-ranef(m, condVar = TRUE)
  qq <- attr(ranef(m, condVar = TRUE)[[1]], "postVar")
  rand.interc<-randoms[[1]] # assuming only one batch of ranefs
  df<-data.frame(Intercepts=rand.interc[,1],
                 sd.interc=2*sqrt(qq[,,1:length(qq)]),
                 lev.names=rownames(rand.interc))
  df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])
  p <- ggplot(df,aes(lev.names,Intercepts))
  p <- p + geom_hline(yintercept=0) +geom_errorbar(aes(ymin=Intercepts-sd.interc, ymax=Intercepts+sd.interc), width=0,color="black") + geom_point(aes(size=1)) 
  p <- p + guides(size=FALSE,shape=FALSE) #+ scale_shape_manual(values=c(1,1,1,16,16,16))
  p <- p + theme_bw() + xlab("") + ylab(axis.label)
  p <- p + theme(axis.text.x=element_text(size=rel(1.2)),
                 axis.title.x=element_text(size=rel(1.3)),
                 axis.text.y=element_text(size=rel(1.2)),
                 panel.grid.minor=element_blank(),
                 panel.grid.major.x=element_blank())
  p <- p+ coord_flip()
  return(p)
}


glmnet.nonpar.boot <- function(x, y, alpha, n.resamples, n.folds, criterion) {
  # run a nonparametric bootstrap on a glmnet regression using explanatory data x and response variable y, bootstrapping n.resamples times, and selecting an alpha value each time using criterion (e.g. "lambda.1se"). alpha is the a parameter for glmnet which selects ridge regression (a=0), LASSO (a=1) or in between.
    coef_samples = matrix(0, n.resamples, ncol(x)+1) # store coefs for each bootstrap iteration -- note have to also include a column for the intercept apparently 
    for (i in 1:n.resamples) {
      z = sample(1:nrow(x), size=nrow(x), replace=TRUE)
      ytemp  = y[z]
      xtemp = x[z,]
      fit1 = cv.glmnet(xtemp, ytemp, alpha=alpha, nfolds=n.folds)
      coef_samples[i,] = matrix(coef(fit1, s=criterion))
    }
    return(coef_samples)
}

get.CIs <- function(coef_samples) {
  z = apply(coef_samples, 2, quantile, probs=c(0.025, 0.5, 0.975))
  return(z)
}



cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}
