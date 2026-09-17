## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,eval=TRUE,fig.align="center",fig.width = 7,fig.height = 5)
old <- options(digits = 2)

## -----------------------------------------------------------------------------
library(vimpclust)
head(HDdata)

## ---- echo = FALSE, eval = TRUE-----------------------------------------------
res <- sparsewkm(X = HDdata[,-14], centers = 2)

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  res <- sparsewkm(X = HDdata[,-14], centers = 2)

## -----------------------------------------------------------------------------
res$W[,1:5]

## -----------------------------------------------------------------------------
plot(res, what="weights.features")

## -----------------------------------------------------------------------------
res$Wm[,1:5]

## -----------------------------------------------------------------------------
plot(res, what="weights.levels")

## -----------------------------------------------------------------------------
plot(res, what="weights.features", Which=c(4,5,11,12))

## -----------------------------------------------------------------------------
plot(res, what="weights.levels", Which=c(11,12))

## -----------------------------------------------------------------------------
plot(res, what="sel.features")
plot(res, what="sel.levels")

## -----------------------------------------------------------------------------
plot(res, what="expl.var")
plot(res, what="w.expl.var")

## ---- message=FALSE-----------------------------------------------------------
library(mclust)
sapply(c(1,5,6), function(x) {adjustedRandIndex(res$cluster[,x],HDdata$HD)})
table(HDdata$HD, res$cluster[,1])
table(HDdata$HD, res$cluster[,5])
table(HDdata$HD, res$cluster[,6])

## -----------------------------------------------------------------------------
info_clust(res, 5, X = HDdata)

## ---- include=FALSE-----------------------------------------------------------
options(old)

