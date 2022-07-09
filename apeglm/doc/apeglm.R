## ----setup, echo=FALSE, results="hide"----------------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=FALSE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  res <- lfcShrink(dds, coef=2, type="apeglm")

## -----------------------------------------------------------------------------
library(airway)
data(airway)
head(assay(airway))

## -----------------------------------------------------------------------------
keep <- head(which(rowSums(assay(airway)) >= 10), 3000)
airway <- airway[keep,]

## -----------------------------------------------------------------------------
library(DESeq2)
dds <- DESeqDataSet(airway, ~cell + dex)
dds$dex <- relevel(dds$dex, "untrt")
dds <- DESeq(dds)
res <- results(dds)

## -----------------------------------------------------------------------------
x <- model.matrix(design(dds), colData(dds))
param <- dispersions(dds)
mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
offset <- matrix(log(sizeFactors(dds)),
                 ncol=ncol(dds),
                 nrow=nrow(dds),byrow=TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  # original interface, also may give Lapack error
#  fit <- apeglm(Y=airway, x=x, log.lik=logLikNB, param=param, coef=ncol(x),
#                threshold=log(2) * 1, mle=mle, offset=offset)

## -----------------------------------------------------------------------------
library(apeglm)
system.time({
  fitR <- apeglm(Y=airway, x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomR")
})
fit <- fitR
names(fit)
str(fit$prior.control)

## -----------------------------------------------------------------------------
system.time({
  fitCR <- apeglm(Y=assay(airway), x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomCR")
})

## -----------------------------------------------------------------------------
system.time({
  fitC <- apeglm(Y=assay(airway), x=x, log.lik=NULL, param=param, coef=ncol(x),
                 threshold=log(2) * 1, mle=mle, offset=offset, method="nbinomC")
})

## -----------------------------------------------------------------------------
class(fit$ranges)
mcols(fit$ranges, use.names=TRUE)

## -----------------------------------------------------------------------------
system.time({
  res.shr <- lfcShrink(dds, coef=5, type="normal")
})
DESeq2.lfc <- res.shr$log2FoldChange
apeglm.lfc <- log2(exp(1)) * fit$map[,5]

## ----apevsdeseq---------------------------------------------------------------
plot(DESeq2.lfc, apeglm.lfc)
abline(0,1)

## ----maplots, fig.width=9, fig.height=3.5-------------------------------------
par(mfrow=c(1,3))
lims <- c(-8,8)
hline <- function() abline(h=c(-4:4 * 2),col=rgb(0,0,0,.2))
xlab <- "mean of normalized counts"
plot(res$baseMean, res$log2FoldChange, log="x",
     ylim=lims, main="MLE", xlab=xlab)
hline()
plot(res$baseMean, DESeq2.lfc, log="x",
     ylim=lims, main="DESeq2", xlab=xlab)
hline()
plot(res$baseMean, apeglm.lfc, log="x",
     ylim=lims, main="apeglm", xlab=xlab)
hline()

## ----pvalcomp, fig.width=8, fig.height=4--------------------------------------
par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(res$pvalue, fit$fsr, col="blue",
     xlab="DESeq2 pvalue", ylab="apeglm local FSR",
     xlim=c(0,1), ylim=c(0,.5))
abline(0,1)
plot(-log10(res$pvalue), -log10(fit$fsr),
     xlab="-log10 DESeq2 pvalue", ylab="-log10 apeglm local FSR",
     col="blue")
abline(0,1)

## ----pvalcomp2----------------------------------------------------------------
plot(res$padj, fit$svalue, col="blue",
     xlab="DESeq2 padj", ylab="apeglm svalue",
     xlim=c(0,.2), ylim=c(0,.02))

## ----maplotthresh-------------------------------------------------------------
s.val <- svalue(fit$thresh)
cols <- ifelse(s.val < .01, "red", "black")
keep <- res$baseMean > 10 # filter for plot
plot(res$baseMean[keep],
     log2(exp(1)) * fit$map[keep,5],
     log="x", col=cols[keep],
     xlab=xlab, ylab="LFC")
abline(h=c(-1,0,1), col=rgb(1,0,0,.5), lwd=2)

## -----------------------------------------------------------------------------
library(emdbook)
n <- 20
f <- factor(rep(1:2,each=n))
mu <- ifelse(res$baseMean > 50, res$baseMean, 50)
set.seed(1)
cts <- matrix(rnbinom(nrow(dds)*2*n, 
                      mu=mu,
                      size=1/dispersions(dds)),
              ncol=2*n)
theta <- runif(nrow(cts),1,1000)
prob <- rnorm(nrow(cts),.5,.05) # close to 0.5
ase.cts <- matrix(rbetabinom(prod(dim(cts)), prob=prob,
                             size=cts, theta=rep(theta,ncol(cts))),
                  nrow=nrow(cts))
idx <- 1:10
idx2 <- which(f == 2)
theta[idx] <- 1000
prob[idx] <- 0.75
# the spiked in genes have an allelic ratio of 0.75
ase.cts[idx,idx2] <- matrix(rbetabinom(length(idx)*length(idx2), prob=prob[idx],
                                       size=cts[idx,idx2], theta=theta[idx]),
                            nrow=length(idx))

## -----------------------------------------------------------------------------
betabinom.log.lik <- function(y, x, beta, param, offset) {
  xbeta <- x %*% beta
  p.hat <- (1+exp(-xbeta))^-1
  dbetabinom(y, prob=p.hat, size=param[-1], theta=param[1], log=TRUE)
}

## -----------------------------------------------------------------------------
theta.hat <- 100 # rough initial estimate of dispersion
x <- model.matrix(~f)
niter <- 3
system.time({
  for (i in 1:niter) {
    param <- cbind(theta.hat, cts)
    fit.mle <- apeglm(Y=ase.cts, x=x, log.lik=NULL, param=param,
                      no.shrink=TRUE, log.link=FALSE, method="betabinCR")
    theta.hat <- bbEstDisp(success=ase.cts, size=cts,
                           x=x, beta=fit.mle$map,
                           minDisp=.01, maxDisp=5000)
  }
})

## ----mlemaplot----------------------------------------------------------------
coef <- 2
xlab <- "mean of normalized counts"
plot(res$baseMean, fit.mle$map[,coef], log="x", xlab=xlab, ylab="log odds")
points(res$baseMean[idx], fit.mle$map[idx,coef], col="dodgerblue", cex=3)

## -----------------------------------------------------------------------------
mle <- cbind(fit.mle$map[,coef], fit.mle$sd[,coef])
param <- cbind(theta.hat, cts)
system.time({
  fit2 <- apeglm(Y=ase.cts, x=x, log.lik=NULL, param=param,
                 coef=coef, mle=mle, threshold=0.5,
                 log.link=FALSE, method="betabinCR")
})

## ----asemaplot, fig.width=8, fig.height=4-------------------------------------
par(mfrow=c(1,2))
ylim <- c(-1,1.5)
s.val <- svalue(fit2$thresh) # small-or-false-sign value
plot(res$baseMean, fit.mle$map[,coef], main="MLE",
     log="x", xlab=xlab, ylab="log odds", ylim=ylim)
points(res$baseMean[idx], fit.mle$map[idx,coef], col="dodgerblue", cex=3)
abline(h=0,col=rgb(1,0,0,.5))
cols <- ifelse(s.val < .01, "red", "black")
plot(res$baseMean, fit2$map[,coef], main="apeglm",
     log="x", xlab=xlab, ylab="log odds", col=cols, ylim=ylim)
points(res$baseMean[idx], fit2$map[idx,coef], col="dodgerblue", cex=3)
abline(h=0,col=rgb(1,0,0,.5))

## -----------------------------------------------------------------------------
logit <- function(x) log(x/(1-x))
logit(.75)
table(abs(logit(prob[s.val < .01])) > .5)

## -----------------------------------------------------------------------------
sessionInfo()

