## ----------[ lodading necessary libraries ]--------------------
library(stats)          #fitting procedure
library(expm)           #Matrix exponential also loads library(Matrix)
library(multicore)      #Parallelisation of code


## ****************************************************************************************
## ----------[ Functions for least squares fitting ]---------------------------------------
## ****************************************************************************************
##' Fits data to aggregate Markov Chain model
##'
##' @title StammFitKstt
##' @param g.dat gene expression data
##' @param t.dat time points of data
##' @param lambda L1 penalty parameter (default = 0.01)
##' @param n.states number of states in the fitting
##' @param fit.as Fitting lin, log2Dat, logDat, log2Al (default = 'lin')
##' @param fix.w Logical variable, fit with fixed "W" matrix only the beta parameter
##' if true (default FALSE)
##' @param w pass the W matrix if
##' @return The function returns fit which is returned from the fitting function.
##' It returns the fitted $w$ matrix and the $/beta$ matrix. It also returns a  obj vector
##' that contains the rss, bic and aic scores for the fit.
##' @author anas ahmad rana
StammFitKstt <- function(g.dat, t.dat, lambda=0.01, n.states=3, fix.w=FALSE, w=NULL, max.fit.iter=4) {
    p <- nrow(g.dat)
    if (fix.w) {
        p <- 1
        x0 <- runif (p * n.states)
        wFit <- w
    } else if (fix.w == FALSE)  {
        x0 <- runif ((n.states - 1) + nrow(g.dat) * n.states )
        wFit <- NULL
    }

    g.dat.l <- asinh(g.dat)
    if ( !is.vector(g.dat) )
        g.nl <- apply(g.dat, 1, sd)
    else
        g.nl <- sd(g.dat)

    fun <- function(x) {
        tmp <- StammPar(x = x, n.states = n.states, p = p, fix.w = fix.w, wFit = wFit)
        wFit <- tmp$w
        betaFit <- tmp$beta
        fit <- StammKstt(wFit, betaFit, t.dat)
        rss <- {asinh(fit$y) - g.dat.l}^2
        penalty <- lambda * sum(betaFit / g.nl)
        ss <- c(rss, penalty)
        obj <- sum(ss)
        obj
    }
    if (is.vector(g.dat)) {
        par.scale <- c(rep(1, n.states - 1), rep(max(g.dat), n.states))
    } else {
        par.scale <- c(rep(1, n.states - 1), rep(apply(g.dat, 1, max), n.states))
    }

    fit.conv <- 10^8
    fit.iter <- 1
    while (fit.conv != 0 & fit.iter < max.fit.iter) {
        res.tmp <- nlminb(x0, fun, lower = 0, upper = max(g.dat), scale = 1 / par.scale,
                  control=list(iter.max=10000, eval.max=7000, rel.tol=10^-14, sing.tol=10^-14))
        if (res.tmp$convergence <= fit.conv){
            fit.conv <- res.tmp$convergence
            res <- res.tmp
        }
        if (res.tmp$convergence > 0)
            writeLines(paste('Not converged, iteration', fit.iter))
        fit.iter  <- fit.iter + 1
    }


    par <- StammPar(g.dat, res$par, n.states, fix.w=fix.w, wFit=w)

    n <- ncol(g.dat)
    if (fix.w) {
        rss <- res$objective - lambda * sum(par$beta / g.nl)
        obj <- rss
        names(obj) <- 'rss'
    } else {
        rss <- res$objective - lambda * sum(par$beta / g.nl)
        obj <- StammBic(rss, n, par$beta)
    }

    return(list(fit = res, w = par$w, beta = par$beta, ms = obj))
}

##' Calculates BIC and AIC
##'
##' BIC and AIC calculated for least squared fit
##' @title StammBic
##' @param rss The residual sum of squares of the fit
##' @param n numer of time points used when fitting
##' @param beta the number of parameters used to calculate Df
##' @return vector of rss bic and aic
##' @author anas ahmad rana
StammBic <- function(rss, n, beta, b.thresh = 10^-4) {
    Df <- sum(beta > b.thresh)
    bicSc <- n * log(rss / (n - 1)) + log(n) * Df
    aicSc <- n * log(rss / (n - 1)) + 2 * Df
    aicScc <- aicSc + (2 * Df * (Df + 1)) / (n - Df - 1)
    ms <- c(rss, bicSc, aicSc, aicScc)
    names(ms) <- c('rss', 'bic', 'aic', 'aicc')
    return(ms)
}

##' Reshapes parameter vecot, x, into beta matrix and w matrix (or only beta matrix)
##'
##' .. content for \details{} ..
##' @title StammPar
##' @param g.dat data matrix used for naming beta
##' @param x parameter vector to reshape into beta and w
##' @param n.states number of states in model
##' @param p no of genes to be fitted (default # rows in data)
##' @param fix.w logical if w is kept fixed or not
##' @return The function returns w only if fix.w=F and it returns beta matrix rearranged from the x vector.
##' @author anas ahmad rana
StammPar <- function(g.dat=NULL, x, n.states, p=nrow(g.dat), fix.w=FALSE, wFit=NULL, fit.cent=FALSE) {
    ## only rearrange
    if (fit.cent) {
        p <- nrow(g.dat)
        betaFit <- matrix(x, p, n.states)
        return(list(w = wFit, beta = betaFit))
    } else if (fix.w && !fit.cent) {
        if (length(x) != n.states )
            stop('No of parameters has to equal no of states when fitting per gene')
        betaFit <- matrix(x, 1, n.states)
        return(list(w = wFit, beta = betaFit))
    } else {
        ## W matrix from x[1:n-1]
        if (n.states == 2) {
            wFit <- matrix(0, n.states, n.states)
            wFit[2, 1] <- x[1]
        } else if (n.states >= 3)  {
            wFit <- matrix(0, n.states, n.states)
            diag(wFit[-1, ]) <- x[1:(n.states - 1)]
        } else {
            wFit <- NULL
        }
        ## Assign other x values to beta
        if (n.states == 1) {
            betaFit <- matrix(x, p, n.states)
        } else {
            lnx <- length(x)
            p <- {lnx - n.states + 1} / n.states
            betaFit <- matrix( x[-c(1:(n.states - 1))], p, n.states)
        }
        if (!is.null(g.dat)) {
            rownames(betaFit) <- rownames(g.dat)
        }
        return(list(w = wFit, beta = betaFit))
    }
}

##' Function that takes as arguments the w_fit matrix, the beta_fit matrix and time
##' points, it then calculates a trajectory
##'
##'
##' @title StammKstt
##' @param wFit W transition matrix default=NULL
##' @param betaFit beta matrix p x T(n)
##' @param t time points
##' @return S unlogged trajectory for the arguments. P state occupation probabilities.
##' @author anas ahmad rana
StammKstt <- function(wFit=NULL, betaFit, t) {
    n <- length(t)
    p <- nrow(betaFit)
    ## create new matrix containing -offdiag on diag
    if (!is.null(wFit)) {
        odiagElements <- diag(wFit[-1, , drop=FALSE ])
        diag(wFit) <- c(-odiagElements, 0)
        k <- nrow(wFit)
        ## initial condition
        p0 <- rep(0, k)
        p0[1] <- 1
        P <- matrix(NA, ncol(betaFit), n)
        for (i in 1:n) {
            P[, i] <- expm::expm(wFit * t[i], method='Ward77') %*% p0
            ## mean gene expression
            S <- betaFit %*% P
        }
    } else {
        S <- rep(betaFit, n)
        P <- 1
    }
    return(list(y = S, prob = P))
}

## ******************************************************************************************
## Clustering and fitting centroids
## ******************************************************************************************

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param g.dat
##' @param t.dat
##' @param k.stt
##' @param m.cl
##' @param n.core
##' @param n.genes
##' @return
##' @author anas ahmad rana
StammCrossValClust.p <- function(g.dat, t.dat, k.stt, m.cl = seq(4, 20, 2),
                                n.core = 20, n.genes=nrow(g.dat)) {
    t.ko <- 2:length(t.dat)

    fit.cl <- vector('list', length(m.cl) * (length(t.ko)))
    fit.gn <- vector('list', length(m.cl) * (length(t.ko)))
    i.a <- 1
    l.name <- NULL
    for (i.m in m.cl) {
        for (i.t in t.ko) {
            print(paste('fitting m =', i.m, 'time deletion ', i.t, '...'))
            g.dat.t <- g.dat[, -i.t]
            g.sd <- apply(g.dat.t, 1, sd)
            g.norm <- g.dat.t / g.sd
            g.km <- kmeans(g.norm, i.m, iter.max=100, nstart=50)
            fit.cl[[i.a]] <- StammFitKstt(g.km$centers, t.dat[-i.t], lambda=0, n.states=k.stt)
            w <- fit.cl[[i.a]]$w
            fit.gn[[i.a]] <- StammFitGns(g.dat=g.dat.t, t.dat=t.dat[-i.t], n.states=k.stt,
                                        w=w, pll=TRUE, n.core=n.core)
            print('... DONE')
            i.a <- i.a + 1
            l.name <- append(l.name, paste('m.', i.m, '_tDel.', i.t, sep=''))
        }
    }
    names(fit.cl) <- l.name
    return(list(fit.cluster=fit.cl, fit.genes=fit.gn))

}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param fit
##' @param t.dat
##' @param g.dat
##' @param t.ko
##' @return
##' @author anas ahmad rana
TDelmse <- function(fit, t.dat, g.dat,  t.ko=2:length(t.dat)) {
    g.sd <- apply(asinh(g.dat), 1, sd)
    i.f <- 1
    mse.vec <- rep(NA, length(t.ko))
    for (i.t in 1:length(t.ko)) {
        beta <- fit[[i.t]]$beta
        w <- fit[[i.t]]$w
        if (is.null(beta))
            stop('wrong variable structure')
        tmp <- StammKstt(w, beta, t.dat[t.ko[i.t]])
        mse.vec[i.t] <- mean(((asinh(tmp$y) - asinh(g.dat[, t.ko[i.t]])) / g.sd)^2)
    }
    mse <- mean(sqrt(mse.vec))
    return(mse)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title StammLoocv.p
##' @param g.dat data to be used when fitting
##' @param t.dat
##' @param k.stt
##' @param m.cl
##' @param n.core
##' @param n.genes
##' @param t.ko
##' @return
##' @author anas ahmad rana
StammLoocv.p <- function(g.dat, t.dat, k.stt, m.cl, n.core=20, n.genes=nrow(g.dat), lambda=0,
                        t.ko=2:length(t.dat)) {
    fit.cl <- vector('list', length(t.ko))
    fit.gn <- vector('list', length(t.ko))
    for (i.t in 1:length(t.ko)) {
        cat(paste('fitting m =', m.cl, 'k =', k.stt, 'time deletion ', t.ko[i.t]))
        g.dat.t <- g.dat[, -t.ko[i.t]]
        g.sd <- apply(g.dat.t, 1, sd)
        g.norm <- g.dat.t / g.sd
        g.km <- kmeans(g.norm, m.cl, iter.max=50, nstart=50)
        cat('.')
        fit.cl[[i.t]] <- StammFitKstt(g.km$centers, t.dat[-t.ko[i.t]], lambda=lambda,
                                     n.states=k.stt)
        cat('.')
        w <- fit.cl[[i.t]]$w
        fit.gn[[i.t]] <- StammFitGns(g.dat=g.dat.t, t.dat=t.dat[-t.ko[i.t]], n.states=k.stt,
                                    lambda=lambda, w=w, pll=TRUE, n.core=n.core)
        cat('. DONE\n')
    }
    names(fit.cl) <- paste('tDel_', t.ko, sep='')
    names(fit.gn) <- paste('tDel_', t.ko, sep='')
    mse <- TDelmse(fit.gn, t.dat, g.dat, t.ko)
    return(list(fit.clust=fit.cl, fit.genes=fit.gn, mse=mse))
}

##' Scale data and calculate k-means cluster for all vector elements
##' of n.cl
##'
##' .. content for \details{} ..
##' @title StammCluster
##' @param g.dat data to be clustered
##' @param n.cl vector of arbitrary length
##' @return
##' @author anas ahmad rana
StammCluster <- function(g.dat, m.cl) {
    ## Scale input data to unit variance
    g.sd <- apply(g.dat, 1, sd)
    g.norm <- (g.dat) / g.sd
    ## initialise empty lists
    g.cl <- vector('list', length(m.cl))
    g.cent <- vector('list', length(m.cl))
    g.cl.names <- vector('list', length(m.cl))
    ## perform a k-means clustering
    for (i.n in 1:length(m.cl)) {
        i.m  <- m.cl[i.n]
        g.km.cl <- kmeans(g.norm, i.m, iter.max=100, nstart=100)
        g.cl[[i.n]] <- g.km.cl
        g.cent[[i.n]] <- g.km.cl$centers
        g.rep.names <- rep(NA, i.m)
        for (i.clg in 1:nrow(g.km.cl$centers)) {
            d.g <- apply((abs(g.norm - g.km.cl$centers[i.clg, ])), 1, sum)
            tmp <- which(d.g == min(d.g))
            g.rep.names[i.clg] <- rownames(g.norm)[tmp]
        }
        g.cl.names[[i.n]] <- g.rep.names
    }
    names(g.cl) <- paste('m.', m.cl, sep='')
    names(g.cent) <- paste('m.', m.cl, sep='')
    names(g.cl.names)  <- paste('m.', m.cl, sep='')
    ## return variables k-means centroid, fit and names of
    ## representative genes
    return(list(cent.dat=g.cent, cl.kmean=g.cl, rep.gns=g.cl.names))
}


##' Cluster (kmeans) and fit under timepoint deletions
##'
##' .. content for \details{} ..
##' @title
##' @param g.dat
##' @param t.dat
##' @param m.cl
##' @param t.ko
##' @param k.stt
##' @return
##' @author anas ahmad rana
ParClusterCV <- function(g.dat, t.dat=NULL, m.cl=seq(5, 20, 2), t.ko=2:ncol(g.dat),
                         k.stt=4, l.pen=0) {
    ## some parameters determined from the input
    n.genes <- nrow(g.dat)
    g.norm <- g.dat / apply(g.dat, 1, sd)
    vec.mt <- cbind(rep(1:length(m.cl), length(t.ko)), rep(t.ko, each=length(m.cl)))

    fit <- mclapply(1:nrow(vec.mt), function(x) {
        t.dl <- vec.mt[x, 2]
        ## Clustered after time point deletion using k-means
        g.km <- kmeans(g.norm[, -t.dl], vec.mt[x, 1], iter.max=100, nstart=20)
        g.ct <- g.km$centers
        ## fit cluster centroids with fixed states k
        fit.m <- StammFitKstt(g.ct, t.dat[-t.dl], lambda=0, n.states=k.stt)
        ## Fix w and fit all genes used in clustering
        w <- fit.m$w
        beta <- matrix(0, n.genes, k.stt)
        rownames(beta) <- rownames(g.dat[1:n.genes, ])
        for (i.j in 1:n.genes) {
            fit.tmp <- StammFitKstt(g.dat[i.j, -t.dl], t.dat[-t.dl], lambda=l.pen,
                                   n.states=k.stt, fix.w=TRUE, w=w)
            beta[i.j, ] <- fit.tmp$beta
        }
        return(list(cl=g.km, cl.fit=fit.m, w=w, beta=beta))
    }, mc.preschedule = TRUE)

    names(fit) <- paste('m', vec.mt[, 1], 'td', vec.mt[, 2], sep='.')
    ## output from the function
    return(list(g.norm=g.norm, fit=fit))
}

StammChsKL <- function (g.dat, t.dat, m, k.vec=2:5, pen.vec=seq(0, 0.2, 0.05), pll=TRUE, w=NULL) {
    n.k <- length(k.vec)
    fit.k <- vector('list', n.k)
    beta <-  vector('list', n.k)
    w <- vector('list', n.k)
    bic <- rep(NA, n.k)
    aicc <- rep(NA, n.k)
    l.min <- rep(NA, n.k)
    for (i.k in 1:n.k) {
        fit.tmp <- FitClGns(g.dat, t.dat, m=m, l.pen=pen.vec, k.stt=k.vec[i.k], pll=pll)
        fit.k[[i.k]] <- list(fit=fit.tmp$fit.g[[which(fit.tmp$bic == min(fit.tmp$bic))]],
                             bic=fit.tmp$bic, aicc=fit.tmp$aicc, all.fit=fit.tmp)
        beta[[i.k]] <- fit.k[[i.k]]$fit$beta
        w[[i.k]] <- fit.k[[i.k]]$fit$w
        bic[[i.k]] <- min(fit.tmp$bic)
        aicc[[i.k]] <- min(fit.tmp$aicc)
        l.min[i.k] <- pen.vec[which(fit.tmp$bic == min(fit.tmp$bic))]
    }
    return(list(k.fit=fit.k, beta=beta, w=w, bic=bic, l.min=l.min))
}

FitClGns <- function(g.dat, t.dat, l.pen=0, k.stt, m, pll=FALSE, w=NULL, n.core=20) {
    if (is.null(w)) {
        ## Normalise data to be univariate and fit clusters
        g.norm <- g.dat / apply(g.dat, 1, sd)
        g.km <- kmeans(g.norm, m, iter.max=1000, nstart=1000)
        g.ct <- g.km$centers
        ## fit cluster centroids with fixed states k
        fit.m <- StammFitKstt(g.ct, t.dat, lambda=0, n.states=k.stt)
        w <- fit.m$w
    }
    ## Fit all genes in g.dat with w from clustering and single penalty
    if (length(l.pen)==1) {
        fit.g <- StammFitGns(g.dat=g.dat, t.dat=t.dat, lambda=l.pen, n.states=k.stt, w=w, pll=pll)
        bic <- fit.g$ms[2]
        aicc <- fit.g$ms[4]
    } else if (is.vector(l.pen)) {
        fit.g <- vector('list', length(l.pen))
        bic <- rep(NA, length(l.pen))
        aicc <- rep(NA, length(l.pen))
        for (i.l in 1:length(l.pen)) {
            fit.g[[i.l]] <- StammFitGns(g.dat=g.dat, t.dat=t.dat, lambda=l.pen[i.l],
                                       n.states=k.stt, w=w, pll=pll, n.core=n.core)
            bic[i.l] <- fit.g[[i.l]]$ms[2]
            aicc[i.l] <- fit.g[[i.l]]$ms[4]

        }
    }
    return(list(bic=bic, aicc=aicc, fit.g=fit.g, fit.m=fit.m, cl.km=g.km))
}

##' Fits betas for a list of genes
##'
##'
##' @title rust.fit.gnlst
##' @param g.names
##' @param g.dat data matrix for fitting
##' @param t.dat time vector of data points
##' @param lambda L1 penalty parameter
##' @param n.states number of states to be fitted
##' @param w transition matrix
##' @param pll logical run in parallel if set to true
##' @return
##' @author anas ahmad rana
StammFitGns <- function(g.names=NULL, g.dat, t.dat, lambda=0, n.states, w, pll=FALSE, n.core=20) {
    if (!is.null(g.names)) {
        g.name.idx <- which(rownames(g.dat)==g.names)
    } else {
        g.name.idx <- 1:nrow(g.dat)
    }

    if (pll) {
        fit.g <- mclapply(g.name.idx, function(x)
                          cl.fit = StammFitKstt(g.dat=g.dat[x, ], t.dat=t.dat, lambda=lambda,
                              n.states=n.states, w=w, fix.w=TRUE),
                          mc.cores = n.core)
    } else {
        fit.g <- lapply(g.name.idx, function(x)
                        cl.fit = StammFitKstt(g.dat=g.dat[x, ], t.dat=t.dat, lambda=lambda,
                            n.states=n.states, w=w, fix.w=TRUE))
    }

    ## take out important information from the fitting data
    betas <- NULL
    rss.v <- NULL
    for (il in 1:length(fit.g)) {
        betas <- rbind(betas, fit.g[[il]]$beta)
        rss.v <- c(rss.v, fit.g[[il]]$ms)
    }

    rownames(betas) <- rownames(g.dat[g.name.idx, ])
    names(rss.v) <- rownames(g.dat[g.name.idx, ])
    rss <- sum(rss.v)
    n <- length(t.dat)
    ms <- StammBic(rss, n, betas)

    return(list(beta=betas, rss.v=rss.v, w=w, ms=ms))
}

FixedParCV <- function(g.dat, t.dat, lambda=0, k.stt, m, n.core=30, t.ko=2:length(t.dat)) {
    fit.cl <- vector('list', length(t.ko))
    fit.gn <- vector('list', length(t.ko))
    for (i.t in 1:length(t.ko)) {
        cat(paste('time deletion ', i.t, '& k,', k.stt, 'm', m, 'lambda', lambda, '...'))
        g.dat.t <- g.dat[, -t.ko[i.t]]
        g.sd <- apply(g.dat.t, 1, sd)
        g.norm <- g.dat.t / g.sd
        g.km <- kmeans(g.norm, m, iter.max=100, nstart=50)
        fit.cl[[i.t]] <- StammFitKstt(g.km$centers, t.dat[-t.ko[i.t]], lambda=0, n.states=k.stt)
        w <- fit.cl[[i.t]]$w
        fit.gn[[i.t]] <- StammFitGns(g.dat=g.dat.t, t.dat=t.dat[-t.ko[i.t]], n.states=k.stt,
                                    w=w, pll=TRUE, n.core=n.core)
        print('... DONE')
    }
    return(list(fit.cluster=fit.cl, fit.gene=fit.gn))
}
