##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param g.dat 
##' @param m.v 
##' @param iter.max 
##' @param nstart 
##' @return 
##' @author anas ahmad rana
StammKmeans <- function(g.dat, m.v=2:30, iter.max=1000, nstart=500) {
    ## Transform data and standardize for kmeans
    g.sdA <- apply(asinh(g.dat), 1, sd)
    g.normA <- asinh(g.dat) / g.sdA

    
    within.ss <- rep(NA, length(m.v))
    for (i.m in 1:length(m.v)) {
        g.km <- kmeans(g.normA, m.v[i.m], iter.max=iter.max, nstart=nstart)
        within.ss[i.m] <- g.km$tot.withinss
    }

    J <- within.ss[-1] / within.ss[1:28]

    kmeans.plot <- StammKmeans.plot(1-J, m.v[-1], x.lab='No. of clusters, m',
                                    y.lab=expression(paste(Delta, "J")), p.title=" ")

    return(list(Delta.J=(1-J), m.vec=m.v[-1], kmeans.plot=kmeans.plot))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param g.dat 
##' @param k.states 
##' @param hat.m 
##' @param n.core 
##' @return 
##' @author anas ahmad rana
StammMSEcv.K <- function(g.dat, t.dat, k.states=1:5, hat.m=13, n.core=50, l.pen=0, return.all=FALSE){
    fit.cv <- vector('list', length(k.states))
    mse <- rep(NA, length(k.states))
    for (i.k in 1:length(k.states)) {
        fit.cv[[i.k]] <- StammLoocv.p(g.dat, t.dat, k.stt=k.states[i.k], m.cl=hat.m, n.core=n.core, lambda=l.pen)
        mse[i.k] <- fit.cv[[i.k]]$mse
    }
    names(fit.cv) <- paste("K=", k.states, sep="")
    names(mse) <- paste("K=", k.states, sep="")
    
    hat.K <- k.states[which(mse == min(mse))]

    if (hat.K != max(k.states)){
        fit <- FitClGns(g.dat, t.dat, l.pen=l.pen, k.stt=hat.K,
                        m=hat.m, pll=TRUE, n.core=n.core)
    } else {
        print("MSEcv has min at max K tried. Consider trying larger K!")
    }
    
    if (return.all)
        return(list(fit.cv=fit.cv, mse.cv=mse, fit=fit, hat.K=hat.K))
    else
        return(list(mse.cv=mse, fit=fit, hat.K=hat.K))
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param m.cl 
##' @param k.states 
##' @return 
##' @author anas ahmad rana
StammMSEcv.m <- function(g.dat, t.dat, m.cl=12:20, k.states=1:5, n.core=50, l.pen=0) {
fit.k <- vector('list', length(k.states))
for (i.k in 1:length(k.states)){
    fit <- vector('list', length(m.cl))
    for (i.m in 1:length(m.cl)) {
        fit[[i.m]] <- FitClGns(g.dat, t.dat, l.pen=l.pen, k.stt=k.states[i.k],
                               m=m.cl[i.m], pll=TRUE, n.core=n.core)
    }
    cat(paste('Fitting k =', k.states[i.k], '... Done\n'))
    names(fit) <- paste("m=", m.cl, sep="")
    fit.k[[i.k]] <- fit
}
names(fit.k) <- paste("K=", k.states, sep="")
return(list(fit=fit.k))
}

