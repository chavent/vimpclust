## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,eval=TRUE,fig.align="center",fig.width = 7,fig.height = 5)
old <- options(digits = 2)

## -----------------------------------------------------------------------------
library(vimpclust)
data(DataMice)
DataMice[1:10, 1:4]

## -----------------------------------------------------------------------------
summary(DataMice$Class.mouse)

## ----  message=FALSE, warning=FALSE-------------------------------------------
names(DataMice)[1:20]
index <- unlist(strsplit(names(DataMice)[1:(dim(DataMice)[2]-5)], split="_"))
index <- as.numeric(index[(1:length(index)%%4==2)])

## ----  message=FALSE, warning=FALSE, echo = T---------------------------------
set.seed(1407)
res.mice <- groupsparsewkm(X = DataMice[,(1:length(index))], centers = 8, index = index, verbose = 1)

## ----  eval = F---------------------------------------------------------------
#  set.seed(1407)
#  res.mice <- groupsparsewkm(X = DataMice[,(1:length(index))], centers = 8, index = index, verbose = 1)

## -----------------------------------------------------------------------------
res.mice$Wg[1:20,1:5]

## -----------------------------------------------------------------------------
plot(res.mice, what="weights.groups")

## -----------------------------------------------------------------------------
res.mice$W[1:20,1:5]

## ----  message=FALSE, warning=FALSE-------------------------------------------
plot(res.mice, what = "weights.features")

## ----  message=FALSE, warning=FALSE-------------------------------------------
plot(res.mice, what = "weights.groups", Which=c(1,2,30))

## ----  message=FALSE, warning=FALSE-------------------------------------------
plot(res.mice, what = "weights.features", Which=c(1,2,30))

## -----------------------------------------------------------------------------
plot(res.mice, what="sel.groups")

## -----------------------------------------------------------------------------
plot(res.mice, what="sel.features")

## -----------------------------------------------------------------------------
plot(res.mice, what="expl.var")

## -----------------------------------------------------------------------------
plot(res.mice, what="w.expl.var")

## ---- message=FALSE-----------------------------------------------------------
sapply(1:length(res.mice$lambda), function(x) {mclust::adjustedRandIndex(res.mice$cluster[,x],DataMice$Class.mouse)})
sapply(1:length(res.mice$lambda), function(x) {mclust::adjustedRandIndex(res.mice$cluster[,x],DataMice$Genotype)})
sapply(1:length(res.mice$lambda), function(x) {mclust::adjustedRandIndex(res.mice$cluster[,x],DataMice$Treatment)})
sapply(1:length(res.mice$lambda), function(x) {mclust::adjustedRandIndex(res.mice$cluster[,x],DataMice$Behaviour)})

## ---- include=FALSE-----------------------------------------------------------
options(old)

