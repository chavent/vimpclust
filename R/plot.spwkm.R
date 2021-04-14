#' @title  Plots from a "spwkm" object 
#' @export
#'
#' @description Produces several graphics to help interpreting a \code{spwkm} object.
#'
#' @param x An object of class \code{spwkm}.
#' @param what A character string indicating which element of \code{x} to be plotted. See section
#' "Details" below for further information.  
#' @param Which A numerical vector indexing the groups or the variables to be displayed. See section
#' "Details" below for further information.   
#' @param xtitle The title of the x-axis.
#' @param ytitle The title of the y-axis.
#' @param title The title of the graphic.
#' @param showlegend A boolean. If \code{showlegend=NULL} (default value), the legend is displayed.
#' @param legendtitle The title of the legend.
#' @param ... Further arguments to the \code{plot} function.
#' 
#' @return \item{p}{an object of class \code{ggplot}.}
#'
#' @details
#' The \code{plot} function allows to represent the regularization paths for a grid of values of \code{lambda}, as well as several quality criteria associated to the 
#' clustering. 
#' 
#' For both \code{groupsparsewkm} and \code{sparsewkm} functions, the following options are available:
#' 
#' If \code{what=weights.features}, the regularization paths for the weights associated to the variables are displayed. 
#' 
#' If \code{what=sel.features}, the graph represents the number of selected variables for each value of the regularization parameter \code{lambda}. In the case of 
#' sparse weighted k-means for mixed data, categorical variables are represented with dotted lines so that one easily identifies them. 
#' 
#' If \code{what=expl.var}, the explained variance (computed as the contribution of the between-class variance to the global variance) is displayed. This criterion is 
#' computed for all variables in the data set, without taking into account the weights of the group or of the variables. 
#' 
#' If \code{what=w.expl.var}, the explained weighted variance is computed. The difference with the criterion above is that the weights of the variables are
#' taken into account in the computation. This leads to a criterion which, for large regularization parameters \code{lambda}, may be computed on one variable only, if 
#' its weight becomes equal to 1 and all the others are discarded. 
#' 
#' If \code{what=pen.crit}, the graph displays the evolution of the penalized criterion, maximized by the algorithm. This criterion writes as 
#' the between-class weighted sum-of-squares, penalized by a group L1-norm. For more details on the mathematical expressions, one may refer to Chavel et al. (2020). 
#' 
#' For the outcome of the \code{groupsparsewkm} function trained on numerical data only, two more options are available:
#' 
#' If \code{what=weights.groups}, the regularization paths for the weights associated to the groups of variables are displayed.
#' 
#' If \code{what=sel.groups}, the graph represents the number of selected groups for each value of the regularization parameter \code{lambda}.
#' 
#' For the outcome of the \code{sparsewkm} function trained on mixed data, two more options are also available:
#' 
#' If \code{what=weights.levels}, the regularization paths for the weights associated to the levels of the categorical variables are displayed. 
#' 
#' If \code{what=sel.levels}, the graph represents the number of selected levels associated to the categorical variables plus the number of selected 
#' numerical variables, for each value of the regularization parameter \code{lambda}.
#' 
#' If the number of groups in \code{groupsparsewkm} or if the number of features in \code{sparsewkm} are too large to have easily interpretable graphics, one may select 
#' some groups or some variables using the argument \code{Which}. Note that when training \code{sparsewkm} on mixed data, the initial order of the variables is changed:
#' after the processing step, numerical variables are displayed first, and categorical second. The indexing provided in \code{Which} should take this into account (see the
#' Examples section). 
#' 
#' @references M., Chavent, J. Lacaille, A. Mourer, and M. Olteanu (2020). 
#' Sparse k-means for mixed data via group-sparse clustering. To appear in ESANN proceedings.
#' 
#' @seealso \code{\link{sparsewkm}}, \code{\link{groupsparsewkm}} 
#' 
#' @examples
#' # sparse weighted k-means on mixed data
#' \donttest{
#' data(HDdata)
#' out <- sparsewkm(X = HDdata[,-14], centers = 2)
#' plot(out, what = "weights.features")
#' plot(out, what = "weights.levels")
#' plot(out, what = "sel.features")
#' plot(out, what = "sel.levels")
#' plot(out, what = "expl.var")
#' plot(out, what = "w.expl.var")
#' plot(out, what = "pen.crit")
#' # plot the regularization paths for first three variables only 
#' plot(out, what = "weights.features", Which=1:3)
#'  
#' # group sparse weighted k-means on numerical data
#' data(iris)
#' index <- c(1, 2, 1, 2)
#' out <- groupsparsewkm(X = iris[,-5], centers = 3, index = index)
#' plot(out, what = "weights.groups")
#' plot(out, what = "weights.features")
#' plot(out, what = "sel.groups")
#' plot(out, what = "sel.features")
#' plot(out, what = "expl.var")
#' plot(out, what = "w.expl.var")
#' plot(out, what = "pen.crit")
#' # plot the regularization paths for the variables in the first group only
#' plot(out, what = "weights.features", Which=1)
#' }
#' @import ggplot2 Polychrome

plot.spwkm <- function(x, what="weights.features", Which=NULL, xtitle =NULL, ytitle = NULL, 
                           title = NULL, showlegend = NULL, 
                           legendtitle = NULL, ...) 
{
  if (!inherits(x, "spwkm")) 
    stop("Use only with \"spwkm\" objects")
  if (x$type=="NumGroupSparse")
  {
    if (what=="weights.levels")
      stop("The implementation does not allow to perform group-sparse k-means on mixed data yet!")
    if (what=="sel.levels")
      stop("The implementation does not allow to perform group-sparse k-means on mixed data yet!")
    if (what=="weights.groups")
    {
      data.to.plot <- data.frame(weights=c(t(x$Wg)), lambda=rep(x$lambda, length(unique(x$index))), 
                                 group=rep(rownames(x$Wg), each=length(x$lambda)),
                                 color=as.factor(rep(1:length(unique(x$index)), 
                                                           each=length(x$lambda))))
      color.count <- length(unique(x$index))
      if (is.null(Which)==FALSE)
      {data.to.plot <- data.to.plot[data.to.plot$color%in%Which,]} else 
        Which <- unique(x$index)
      get.palette <-  grDevices::colorRampPalette(glasbey.colors(n=32))
      colors.groups <- get.palette(color.count)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Group weights"
      if (is.null(title)) 
        title = "Groups of features - regularization paths"
      if(is.null(showlegend))
        showlegend = T
      if (is.null(legendtitle)) 
        legendtitle = "Group of features"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$weights, group=as.factor(.data$group))) +
        geom_line(aes(color=.data$color))+
        geom_point(aes(color=.data$color))+
        scale_colour_manual("", values=colors.groups[Which], labels=rownames(x$Wg)[Which]) +
        xlab(xtitle)+ylab(ytitle)+
        labs(color = legendtitle, linetype = legendtitle)+ggtitle(title)
      if (showlegend==F)
        p <- p + theme(legend.position="none")
        return(p)
    }
    if (what=="weights.features")
    {
      data.to.plot <- data.frame(weights=c(t(x$W)), lambda=rep(x$lambda, dim(x$W)[1]), 
                                 group=rep(rownames(x$W), each=length(x$lambda)), 
                                 color=as.factor(rep(x$index, 
                                                     each=length(x$lambda))))
      color.count <- length(unique(x$index))
      if (is.null(Which)==FALSE)
      {data.to.plot <- data.to.plot[data.to.plot$color%in%Which,]} else 
        Which <- unique(x$index)
      get.palette <-  grDevices::colorRampPalette(glasbey.colors(n=32))
      colors.groups <- get.palette(color.count)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Features weights"
      if (is.null(title)) 
        title = "Features - regularization paths"
      if(is.null(showlegend))
        showlegend = T
      if (is.null(legendtitle)) 
        legendtitle = "Features"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$weights, group=as.factor(.data$group), 
                                  color=.data$color)) +
        geom_line()+ geom_point()+
        scale_colour_manual("", values=colors.groups[unique(data.to.plot$color)], 
                            labels= paste("Group ", unique(data.to.plot$color), sep=""))+
        xlab(xtitle)+ ylab(ytitle)+ggtitle(title)
      if (showlegend==F)
        p <- p + theme(legend.position="none")
      return(p)
    }
    if (what=="sel.groups")
    {
      data.to.plot <- data.frame(sel.groups=x$sel.groups, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Number of selected groups"
      if (is.null(title)) 
        title = "Selected groups path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$sel.groups)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="sel.features")
    {
      data.to.plot <- data.frame(sel.features=x$sel.feat, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Number of selected features"
      if (is.null(title)) 
        title = "Selected features path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$sel.features)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="expl.var")
    {
      data.to.plot <- data.frame(bss=apply(x$bss.per.feature,2,sum)/dim(x$W)[1], lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Explained variance"
      if (is.null(title)) 
        title = "Explained variance path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$bss)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="pen.crit")
    {
      penalized.criterion <- sapply(1:length(x$lambda), function(i){sum(x$W[,i]*x$bss.per.feature[,i])-
          x$lambda[i]*sum(sqrt(summary(as.factor(x$index)))*x$Wg[,i])})
      data.to.plot <- data.frame(penalized.criterion, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Penalized weigthed between-class variance"
      if (is.null(title)) 
        title = "Penalized criterion path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$penalized.criterion)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="w.expl.var")
    {
      bss <- sapply(1:length(x$lambda), function(i){sum(x$W[,i]*x$bss.per.feature[,i])})/apply(x$W,2,sum)
      data.to.plot <- data.frame(bss, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Explained weighted variance"
      if (is.null(title)) 
        title = "Explained weighted-variance path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$bss)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
  }
  
  if (x$type=="MixedSparse")
  {
    if (what=="weights.groups")
      stop("The implementation does not allow to perform group-sparse k-means on mixed data yet!")
    if (what=="sel.groups")
      stop("The implementation does not allow to perform group-sparse k-means on mixed data yet!")
    if (what=="weights.features")
    {
      data.to.plot <- data.frame(weights=c(t(x$W)), lambda=rep(x$lambda, length(unique(x$index))), 
                                 feature=rep(rownames(x$W), each=length(x$lambda)),
                                 color=as.factor(rep(1:length(unique(x$index)), 
                                                     each=length(x$lambda))),
                                 line.type=as.factor(rep(as.numeric(summary(as.factor(x$index))>1)*2+1,
                                                     each=length(x$lambda))))
        
      color.count <- length(unique(x$index))
      line.groups <- as.numeric(summary(as.factor(x$index))>1)*2+1
      if (is.null(Which)==FALSE)
        {data.to.plot <- data.to.plot[data.to.plot$color%in%Which,]} else 
        Which <- unique(x$index)
      get.palette <-  grDevices::colorRampPalette(glasbey.colors(n=32))
      colors.groups <- get.palette(color.count)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Features weights"
      if (is.null(title)) 
        title = "Features - regularization paths"
      if(is.null(showlegend))
        showlegend = T
      if (is.null(legendtitle)) 
        legendtitle = ""
      p <- ggplot(data.to.plot, aes(x=.data$lambda, y=.data$weights, 
                               linetype=interaction(.data$color, .data$line.type),
                               col=interaction(.data$line.type, .data$color))) +
        geom_line() +
        geom_point() +
        scale_colour_manual("", values=colors.groups[Which], labels=rownames(x$W)[Which]) +
        scale_linetype_manual("", values=line.groups[Which], labels=rownames(x$W)[Which])+
        xlab(xtitle)+ylab(ytitle)+
        labs(title=legendtitle)+ggtitle(title)
      if (showlegend==F)
        p <- p + theme(legend.position="none")
      return(p)
       
    }
    if (what=="weights.levels")
    {
      data.to.plot <- data.frame(weights=c(t(x$Wm)), lambda=rep(x$lambda, dim(x$Wm)[1]), 
                                 group=rep(rownames(x$Wm), each=length(x$lambda)), 
                                 feature=as.factor(rep(x$index, 
                                                     each=length(x$lambda))),
                                 color=as.factor(rep(1:length(x$index), 
                                                     each=length(x$lambda))))
      if (is.null(Which)==FALSE)
      {data.to.plot <- data.to.plot[data.to.plot$feature%in%Which,]} else
        Which=unique(x$index)
      data.to.plot <- data.to.plot[data.to.plot$feature%in%which(summary(as.factor(x$index))>1),]
      color.count <- length(x$index)
      get.palette <-  grDevices::colorRampPalette(glasbey.colors(n=32))
      colors.groups <- get.palette(color.count)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Levels weights"
      if (is.null(title)) 
        title = "Levels - regularization paths"
      if(is.null(showlegend))
        showlegend = T
      if (is.null(legendtitle)) 
        legendtitle = ""
      p <- ggplot(data.to.plot, aes(x=.data$lambda, y=.data$weights,  group=as.factor(.data$group),
                                    col=.data$color)) +
        geom_line(linetype="dotted") +
        geom_point() +
        scale_colour_manual("", values=colors.groups[which(summary(data.to.plot$color)>0)], 
                            labels= unique(data.to.plot$group))+
        xlab(xtitle)+ylab(ytitle)+
        labs(title=legendtitle)+ggtitle(title)
      if (showlegend==F)
        p <- p + theme(legend.position="none")
      return(p)
    }
    if (what=="sel.features")
    {
      data.to.plot <- data.frame(sel.groups=x$sel.init.feat, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Number of selected features"
      if (is.null(title)) 
        title = "Selected features path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$sel.groups)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="sel.levels")
    {
      data.to.plot <- data.frame(sel.features=x$sel.trans.feat, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Number of selected levels"
      if (is.null(title)) 
        title = "Selected levels path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$sel.features)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="expl.var")
    {
      data.to.plot <- data.frame(bss=apply(x$bss.per.feature,2,sum)/dim(x$Wm)[1], lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Explained variance"
      if (is.null(title)) 
        title = "Explained variance path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$bss)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="pen.crit")
    {
      penalized.criterion <- sapply(1:length(x$lambda), function(i){sum(x$Wm[,i]*x$bss.per.feature[,i])-
          x$lambda[i]*sum(sqrt(summary(as.factor(x$index)))*x$W[,i])})
      data.to.plot <- data.frame(penalized.criterion, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Penalized weigthed between-class variance"
      if (is.null(title)) 
        title = "Penalized criterion path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$penalized.criterion)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
    if (what=="w.expl.var")
    {
      bss <- sapply(1:length(x$lambda), function(i){sum(x$Wm[,i]*x$bss.per.feature[,i])})/apply(x$Wm,2,sum)
      data.to.plot <- data.frame(bss, lambda=x$lambda)
      if(is.null(xtitle))
        xtitle = "Lambda"
      if(is.null(ytitle))
        ytitle = "Explained weighted variance"
      if (is.null(title)) 
        title = "Explained weighted-variance path"
      p<-ggplot(data.to.plot, aes(x=.data$lambda, y=.data$bss)) +
        geom_line()+geom_point()+xlab(xtitle)+ylab(ytitle)+ggtitle(title)
      return(p)
    }
  }
}

