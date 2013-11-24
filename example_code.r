library(stamm)

data(gdat)
data(tdat)
data(sim_data)


fit.km <- StammKmeans(g.dat, m.v=2:30, iter.max=1000, nstart=500)
print(fit.km$kmeans.plot)

StammMSEcv.K(g.dat, t.dat, k.states=1:5, hat.m=13, n.core=50, l.pen=0, return.all=FALSE)
StammMSEcv.m(g.dat, t.dat, m.cl=12:20, k.states=1:5, n.core=50, l.pen=0)

StammMSEcvK.plot(mse, k.states)
StammStab.plot(fit.k, m.v, m.init=2)
StammStabPen.plot(fit, fit.m, t.dat)
