library(stamm)

data(gdat)
data(tdat)
data(sim_data)


fit.km <- StammKmeans(g.dat, m.v=2:30)
print(fit.km$kmeans.plot)

k.stt <- 1:5
m.cl <- 12:20

fit.K <- StammMSEcv.K(g.dat, t.dat, k.states=k.stt, hat.m=13, n.core=50, l.pen=0, return.all=TRUE)
fit.m <- StammMSEcv.m(g.dat, t.dat, m.cl=m.cl, k.states=k.stt, n.core=50, l.pen=0)

StammMSEcvK.plot(fit.K$mse.cv, k.states=1:5)
StammStab.plot(fit.m$fit, m.v=m.cl, m.init=2, k.states=k.stt)
i.k <- which(k.stt==fit.K$hat.K)
StammStabPen.plot(fit.K$fit, fit.K$fit.cv[[i.k]], t.dat)
