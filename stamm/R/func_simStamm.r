library(msm)
library(tmvtnorm)

## ----------[ Simulations on single cells]--------------------
##' Forward simulation of model using single cells
##'
##' @title StammSim
##' @param n.cells number of cells used in simulation (default 200)
##' @param n.genes number of genes simulated (default 9)
##' @param tau average jump.time per state
##' @param n.states number of states simulated (default 2)
##' @param beta.vec predefined beta values for simulation (default randomly generated)
##' @param dt timestep used throughout simulation
##' @param end.time final timepoint (in "real" time)
##' @param av.noise standard deviation of Gaussian noise added to population average (default 0.01)
##' @param stt.noise vector of standard deviations for Gaussian noise per state (default 0)
StammSim <- function(n.cells = 200, n.genes = 9, tau = c(3.5, 5, 14.5), n.states = 2,
                    beta.vec = NULL, dt = 0.01, end.time = 30, av.noise = 0.01,
                    stt.noise = rep(0, n.states), jump.dist = 'exp', sd.par = NULL) {

    ## Checking if some of the arguments passed are the right format and adjusting if possible
    if (!is.vector(beta.vec) & !is.null(beta.vec)) {
        stop('beta must be passed as vector')
    }
    if (length(stt.noise) != n.states & is.vector(stt.noise)) {
        stt.noise <- rep(stt.noise, n.states)
    }
    ## get jump time from the state occupation time use function JumpTime
    if (jump.dist != 'exp' & is.null(sd.par))
        sd.par  <- end.time / n.states
    jump.time <- JumpTime(tau, n.states, n.cells, jump.dist = jump.dist, sd.par = sd.par)
    ## Checks if beta values are passed in the argument (as vector)
    ## Assigns randomvalues to beta if not passed as argument
    if (is.null(beta.vec)) {
        beta.vec <- rnorm(n.genes * n.states, sd = 5)
        beta.vec[beta.vec < 0]  <- 10^(-40)
    }
    ## Reshape betavector as matrix
    betaVals <- matrix(beta.vec, n.genes, n.states)
    ## Initialise results matrix with zeros
    gSim <- matrix(rep(0, n.genes * end.time / dt), n.genes, end.time / dt)
    nPoints <- end.time / dt
    for (iCells in 1:n.cells) {
        gSimC <- NULL
        for (jStates in 1:(n.states - 1)) {
            nSentries <- ceiling(jump.time[jStates, iCells] / dt)
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (nPoints - ncol(gSimC))
        if (nSentries > 0) {
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        gSim <- gSim + gSimC[, 1:nPoints] / n.cells
    }
    ## Add gaussian noise to the data
    datasim <- asinh(gSim) + matrix(rnorm(length(gSim), sd = av.noise), dim(gSim))
    dataSim <- sinh(datasim)
    ## Return values are full simulated data all time points, beta and if t given gSim
    ## with t-pts return all parameters used for this simulation
    return(list(gsim = gSim, beta = betaVals, dataSim = dataSim, n.cells = n.cells, n.gns = n.genes,
                tau = tau, dt = dt, n.stt = n.states, ns.av = av.noise))
}

##'  Function that add noise to normal data
##'
##' @title AddNoise
##' @param sim
##' @param ns.sd
##' @param ns.type
##' @param gData simulated "expression" data
##' @return datasim
##' @author anas ahmad rana
AddNoise <- function(sim = NULL, ns.sd, gData = NULL) {
    if (is.null(gData))
        gData <- sim$gsim
    datsim <- asinh(gData) + matrix(rnorm(length(gData), sd = ns.sd), dim(gData))
    datasim <- sinh(datsim)
    return(datasim)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title JumpTime
##' @param tau average transition times
##' @param n.states number of states simulated
##' @param n.cells number of cells simulated
##' @param jump.dist Distribution of jump time (default = 'exp')
##' @param sd.par second parameter to use for non exponential distributions
##' @return jump.time
##' @author anas ahmad rana
JumpTime  <- function(tau, n.states, n.cells, jump.dist = 'exp', sd.par = NULL) {
    jump.time <- NULL
    if (jump.dist == 'exp') {
        for (i in 1:(n.states - 1)) {
            jump.time <- rbind(jump.time, rexp(n.cells, 1 / tau[i]))
        }
    } else if (jump.dist == 't.norm') {
        for (i in 1:(n.states -1)) {
            jump.time  <- rbind(jump.time, abs(rtnorm(n = n.cells, mean = tau[i], sd = sd.par)))
        }
    } else if (jump.dist == 't.student') {
        jump.time <- matrix(NA, length(tau), n.cells)
        for (i in 1:length(tau)) {
            jump.time[i, ] <- rtmvt(n=n.cells, mean=tau[i], df=sd.par, lower=0)
        }
    }
    return(jump.time)
}

CellFate  <- function(tau, n.cells) {
    cell.fate <- matrix(NA, length(tau), n.cells)
    for (i in 1:length(tau)) {
        cell.fate[i, ] <- rexp(n.cells, 1 / tau[i])
    }
    return(cell.fate)
}

StammSim.rd <- function(n.cells=200, n.genes=9, tau=c(5, 8, 15), n.states=3, beta.vec=NULL, dt=0.01,
                       end.time=30, av.noise=0.01, stt.noise=rep(0, n.states), jump.dist='exp', sd.par=NULL,
                       p.dup=0.1, p.dead=0.1, p.n=0.9) {

    ## Checking if some of the arguments passed are the right format and adjusting if possible
    if (!is.vector(beta.vec) & !is.null(beta.vec)) {
        stop('beta must be passed as vector')
    }
    if (length(stt.noise) != n.states & is.vector(stt.noise)) {
        stt.noise <- rep(stt.noise, n.states)
    }
    ## get jump time from the state occupation time use function JumpTime
    if (jump.dist != 'exp' & is.null(sd.par))
        sd.par  <- end.time / n.states
    jump.time <- JumpTime(tau, n.states, n.cells, jump.dist = jump.dist, sd.par = sd.par)
    ## Checks if beta values are passed in the argument (as vector)
    ## Assigns randomvalues to beta if not passed as argument
    if (is.null(beta.vec)) {
        beta.vec <- rnorm(n.genes * n.states, sd = 5)
        beta.vec[beta.vec < 0]  <- 10^(-40)
    }
    ## Reshape betavector as matrix
    betaVals <- matrix(beta.vec, n.genes, n.states)
    ## Initialise results matrix with zeros
    n.point <- end.time / dt

    gSim <- array(NA, dim=c(n.cells, n.genes, n.point))
    for (iCells in 1:n.cells) {
        gSimC <- NULL
        for (jStates in 1:(n.states - 1)) {
            nSentries <- ceiling(jump.time[jStates, iCells] / dt)
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (n.point - ncol(gSimC))
        if (nSentries > 0) {
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        ## gSim <- gSim + gSimC[, 1:n.point] / n.cells
        gSim[iCells, , ] <- gSimC[, 1:n.point]
    }

    branch.prob <- c(p.dup, p.n, p.dead)
    cell.branch <- sample(c(0,1,2), n.cells,replace=T, prob=branch.prob)

    i.dup <- which(cell.branch==2)
    i.dead <- which(cell.branch==0)
    t.dup <- sample(1:n.point, length(i.dup))
    t.dead <- sample(1:n.point, length(i.dead))

    ## set dead cell to -1 after
    for (i.d in 1:length(i.dead)) {
        gSim[i.dead[i.d], ,-c(1:t.dead[i.d])] <- -1
    }



    jump.dup <- jump.time[, i.dup]
    i.stt.dup <- rep(1, length(i.dup))
    for (i.stt in 1:nrow(jump.dup)) {
        i.update <- t.dup * dt > jump.dup[i.stt, ]
        i.stt.dup[i.update] <- i.stt + 1
    }

    gSim.dup <- array(-1, dim=c(length(i.dup), n.genes, n.point))
    jump.time.dup <- JumpTime(tau, n.states, length(t.dup), jump.dist=jump.dist, sd.par=sd.par)
    for (i.cell in 1:length(i.dup)) {
        gSimC <- NULL
        for (jStates in i.stt.dup[i.cell]:(n.states - 1)) {
            if (jStates < n.states)
                nSentries <- ceiling(jump.time.dup[jStates, i.cell] / dt)
            else
                nSentries <- n.point - t.dup[i.cell]
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (n.point - t.dup[i.cell] + ncol(gSimC))
        if (nSentries > 0) {
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        ## gSim <- gSim + gSimC[, 1:n.point] / n.cells
        gSim.dup[i.cell, , -c(1:t.dup[i.cell])] <- gSimC[, 1:(n.point - t.dup[i.cell])]
    }


    gSim.A <- array(NA, dim=c(n.cells + length(t.dup), n.genes, n.point))
    gSim.A[1:n.cells, , ] <- gSim
    gSim.A[-c(1:n.cells), , ] <- gSim.dup
    tmp.n <- (gSim.A[, 1, ] >= 0)
    n.cells <- apply(tmp.n, 2, sum)

    gSim.A[gSim.A < 0] <- NA

    gSim.mu <- matrix(NA, n.genes, n.point)
    for (i.g in 1:n.genes) {
        gSim.mu[i.g, ] <- apply(gSim.A[, i.g, ], 2, function(x) sum(x, na.rm=TRUE)) / n.cells
    }

    ## Add gaussian noise to the data
    datasim <- asinh(gSim.mu) + matrix(rnorm(length(gSim.mu), sd = av.noise), dim(gSim.mu))
    dataSim <- sinh(datasim)
    ## Return values are full simulated data all time points, beta and if t given gSim
    ## with t-pts return all parameters used for this simulation
    return(list(gsim=gSim.mu, beta=betaVals, dataSim=dataSim, n.cells=n.cells, n.gns=n.genes,
                tau=tau, dt=dt, n.stt=n.states, ns.av=av.noise, cell.branch=cell.branch))
}

StammSim.dd.exmp <- function(n.cells=200, n.genes=9, tau=c(5, 8, 15), n.states=3, beta.vec=NULL,
                            dt=0.01, end.time=30, av.noise=0.01, stt.noise=rep(0, n.states),
                            jump.dist='exp', sd.par=NULL, rate.dup=1/100, rate.dead=1/100) {

    ## Checking if some of the arguments passed are the right format and adjusting if possible
    if (!is.vector(beta.vec) & !is.null(beta.vec)) {
        stop('beta must be passed as vector')
    }
    if (length(stt.noise) != n.states & is.vector(stt.noise)) {
        stt.noise <- rep(stt.noise, n.states)
    }
    ## get jump time from the state occupation time use function JumpTime
    if (jump.dist != 'exp' & is.null(sd.par))
        sd.par  <- end.time / n.states
    jump.time <- JumpTime(tau, n.states, n.cells, jump.dist = jump.dist, sd.par = sd.par)
    ## Checks if beta values are passed in the argument (as vector)
    ## Assigns randomvalues to beta if not passed as argument
    if (is.null(beta.vec)) {
        beta.vec <- rnorm(n.genes * n.states, sd = 5)
        beta.vec[beta.vec < 0]  <- 10^(-40)
    }
    ## Reshape betavector as matrix
    betaVals <- matrix(beta.vec, n.genes, n.states)
    ## Initialise results matrix with zeros
    n.point <- end.time / dt

    gSim <- array(NA, dim=c(n.cells, n.genes, n.point))
    for (iCells in 1:n.cells) {
        gSimC <- NULL
        for (jStates in 1:(n.states - 1)) {
            nSentries <- ceiling(jump.time[jStates, iCells] / dt)
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (n.point - ncol(gSimC))
        if (nSentries > 0) {
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        ## gSim <- gSim + gSimC[, 1:n.point] / n.cells
        gSim[iCells, , ] <- gSimC[, 1:n.point]
    }

    cell.fate <- CellFate(1/c(rate.dup, rate.dead), n.cells)
    i.dup <- which((cell.fate[1,] < end.time) & (cell.fate[1, ] < cell.fate[2,]) )
    i.dead <- which((cell.fate[1, ] > cell.fate[2,]) & cell.fate[2, ] < end.time)
    t.dup <- sample(1:n.point, length(i.dup))
    t.dead <- sample(1:n.point, length(i.dead))

    ## set dead cell to -1 after
    for (i.d in 1:length(i.dead)) {
        gSim[i.dead[i.d], ,-c(1:t.dead[i.d])] <- -1
    }
    ## set duplicated cells to -1 after duplication
    for (i.d in 1:length(i.dup)) {
        gSim[i.dup[i.d], ,-c(1:t.dup[i.d])] <- -1
    }


    jump.dup <- jump.time[, i.dup]
    i.stt.dup <- rep(1, length(i.dup))
    for (i.stt in 1:nrow(jump.dup)) {
        i.update <- t.dup * dt > jump.dup[i.stt, ]
        i.stt.dup[i.update] <- i.stt + 1
    }

    t.dup <- rep(t.dup, 2)
    i.dup <- rep(i.dup, 2)
    i.stt.dup <- rep(i.stt.dup, 2)
    gSim.dup <- array(-1, dim=c(length(i.dup), n.genes, n.point))
    jump.time.dup <- JumpTime(tau, n.states, 2*length(t.dup), jump.dist=jump.dist, sd.par=sd.par)
    for (i.cell in 1:length(i.dup)) {
        gSimC <- NULL
        for (jStates in i.stt.dup[i.cell]:(n.states - 1)) {
            if (jStates < n.states)
                nSentries <- ceiling(jump.time.dup[jStates, i.cell] / dt)
            else
                nSentries <- n.point - t.dup[i.cell]
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (n.point - t.dup[i.cell] + ncol(gSimC))
        if (nSentries > 0) {
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        ## gSim <- gSim + gSimC[, 1:n.point] / n.cells
        gSim.dup[i.cell, , -c(1:t.dup[i.cell])] <- gSimC[, 1:(n.point - t.dup[i.cell])]
    }


    gSim.A <- array(NA, dim=c(n.cells + length(t.dup), n.genes, n.point))
    gSim.A[1:n.cells, , ] <- gSim
    gSim.A[-c(1:n.cells), , ] <- gSim.dup
    tmp.n <- (gSim.A[, 1, ] >= 0)
    n.cell.dup <- apply(tmp.n, 2, sum)

    gSim.A[gSim.A < 0] <- NA

    gSim.mu <- matrix(NA, n.genes, n.point)
    for (i.g in 1:n.genes) {
        gSim.mu[i.g, ] <- apply(gSim.A[, i.g, ], 2, function(x) sum(x, na.rm=TRUE)) / n.cell.dup
    }

    ## Add gaussian noise to the data
    datasim <- asinh(gSim.mu) + matrix(rnorm(length(gSim.mu), sd = av.noise), dim(gSim.mu))
    dataSim <- sinh(datasim)
    ## Return values are full simulated data all time points, beta and if t given gSim
    ## with t-pts return all parameters used for this simulation
    return(list(gsim=gSim.mu, beta=betaVals, dataSim=dataSim, n.cells=n.cell.dup, n.gns=n.genes,
                tau=tau, dt=dt, n.stt=n.states, ns.av=av.noise, cell.fate=cell.fate))
}

StammSim.back <- function(n.cells=200, n.genes=9, tau=c(5, 8, 15), n.states=4, beta.vec=NULL, dt=0.01,
                         end.time=50, av.noise=0.01, stt.noise=rep(0, n.states), jump.dist='exp', sd.par=NULL,
                         back.tau) {

    ## Checking if some of the arguments passed are the right format and adjusting if possible
    if (!is.vector(beta.vec) & !is.null(beta.vec)) {
        stop('beta must be passed as vector')
    }
    if (length(stt.noise) != n.states & is.vector(stt.noise)) {
        stt.noise <- rep(stt.noise, n.states)
    }
    ## get jump time from the state occupation time use function JumpTime
    if (jump.dist != 'exp' & is.null(sd.par))
        sd.par  <- end.time / n.states
    jump.time <- JumpTime(tau, n.states, n.cells, jump.dist = jump.dist, sd.par = sd.par)
    ## Checks if beta values are passed in the argument (as vector)
    ## Assigns randomvalues to beta if not passed as argument
    if (is.null(beta.vec)) {
        beta.vec <- rnorm(n.genes * n.states, sd = 5)
        beta.vec[beta.vec < 0]  <- 10^(-40)
    }
    ## Reshape betavector as matrix
    betaVals <- matrix(beta.vec, n.genes, n.states)
    ## Initialise results matrix with zeros
    n.point <- end.time / dt

    i.stt.fin <- rep(NA, n.cells)
    gSim <- array(NA, dim=c(n.cells, n.genes, n.point))
    for (iCells in 1:n.cells) {
        gSimC <- NULL
        for (jStates in 1:(n.states - 1)) {
            nSentries <- ceiling(jump.time[jStates, iCells] / dt)
            val_tmp <- betaVals[, rep.int(jStates, nSentries)]
            gSimC <- cbind(gSimC,  val_tmp + rnorm(length(val_tmp), sd = stt.noise[jStates]) )
            rm(val_tmp)
        }
        nSentries <- (n.point - ncol(gSimC))
        if (nSentries > 0) {
            i.stt.fin[iCells] <- n.point - ncol(gSimC)
            val_fin <- betaVals[, rep.int(n.states, nSentries)]
            gSimC <- cbind(gSimC, val_fin + rnorm(length(val_fin), sd = stt.noise[n.states]))
            rm(val_fin)
        }
        ## gSim <- gSim + gSimC[, 1:n.point] / n.cells
        gSim[iCells, , ] <- gSimC[, 1:n.point]
    }

    tmp.time <- n.point - i.stt.fin
    i.cellFin <- !is.na(i.stt.fin)
    i.cellFin <- (1:n.cells)[i.cellFin]

    n.rep <- 5
    jump.back <- JumpTime(rep(c(back.tau, tau[n.states-1]), n.rep), 5, n.cells)
    stt.back <- rep(3:4, n.rep)
    for (i.cell in (i.cellFin)) {
        for (i.k in 1:nrow(jump.back)) {
            n.new <- ceiling(jump.back[i.k, i.cell] / dt)
            if (n.new < tmp.time[i.cell]) {
                gSim[i.cell, ,-c(1:i.stt.fin[i.cell])] <- betaVals[, rep.int(stt.back[i.k], length(n.new))]
                tmp.time[i.cell] <- tmp.time[i.cell] - n.new
                i.stt.fin[i.cell] <- i.stt.fin[i.cell] + n.new
            }
        }
    }

    gSim.mu <- matrix(NA, n.genes, n.point)
    for (i.g in 1:n.genes) {
        gSim.mu[i.g, ] <- apply(gSim[, i.g, ], 2, function(x) mean(x, na.rm=TRUE))
    }


    ## Add gaussian noise to the data
    datasim <- asinh(gSim.mu) + matrix(rnorm(length(gSim.mu), sd = av.noise), dim(gSim.mu))
    dataSim <- sinh(datasim)
    ## Return values are full simulated data all time points, beta and if t given gSim
    ## with t-pts return all parameters used for this simulation
    return(list(gsim=gSim.mu, beta=betaVals, dataSim=dataSim, n.cells=n.cells, n.gns=n.genes,
                tau=tau, dt=dt, n.stt=n.states, ns.av=av.noise))
}
