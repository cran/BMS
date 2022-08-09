### R code from vignette source 'bmsmanual.Rnw'

###################################################
### code chunk number 1: bmsmanual.Rnw:16-17
###################################################
options(width=75)


###################################################
### code chunk number 2: bmsmanual.Rnw:64-65
###################################################
data(attitude)


###################################################
### code chunk number 3: bmsmanual.Rnw:68-69
###################################################
library(BMS)


###################################################
### code chunk number 4: bmsmanual.Rnw:72-73
###################################################
att = bms(attitude, mprior = "uniform", g="UIP", user.int=F)


###################################################
### code chunk number 5: bmsmanual.Rnw:79-80
###################################################
coef(att)


###################################################
### code chunk number 6: bmsmanual.Rnw:87-88
###################################################
coef(att, std.coefs=T, order.by.pip=F, include.constant=T) 


###################################################
### code chunk number 7: bmsmanual.Rnw:94-95
###################################################
summary(att) 


###################################################
### code chunk number 8: bmsmanual.Rnw:100-101
###################################################
topmodels.bma(att)[,1:3]


###################################################
### code chunk number 9: bmsmanual.Rnw:106-107
###################################################
image(att)


###################################################
### code chunk number 10: bmsmanual.Rnw:115-116
###################################################
sum(coef(att)[,1])


###################################################
### code chunk number 11: bmsmanual.Rnw:120-121
###################################################
plotModelsize(att)


###################################################
### code chunk number 12: bmsmanual.Rnw:132-133
###################################################
att_fixed = bms(attitude, mprior="fixed", mprior.size=2, user.int=T)


###################################################
### code chunk number 13: bmsmanual.Rnw:139-140
###################################################
att_pip = bms(attitude, mprior="pip", mprior.size=c(.01,.5,.5,.5,.5,.5), user.int=F)


###################################################
### code chunk number 14: bmsmanual.Rnw:146-147
###################################################
plotModelsize(att_fixed) 


###################################################
### code chunk number 15: bmsmanual.Rnw:153-155
###################################################
att_random = bms(attitude, mprior="random", mprior.size=3, user.int=F)
plotModelsize(att_random)


###################################################
### code chunk number 16: bmsmanual.Rnw:161-162 (eval = FALSE)
###################################################
## plotComp(Uniform=att, Fixed=att_fixed, PIP=att_pip, Random=att_random)


###################################################
### code chunk number 17: bmsmanual.Rnw:165-166
###################################################
plotComp(Uniform=att, Fixed=att_fixed, PIP=att_pip, Random=att_random, cex=2)


###################################################
### code chunk number 18: bmsmanual.Rnw:199-201 (eval = FALSE)
###################################################
## data(datafls)
## fls1 = bms(datafls, burn=50000, iter=100000, g="BRIC", mprior="uniform", nmodel=2000, mcmc="bd", user.int=F)


###################################################
### code chunk number 19: bmsmanual.Rnw:203-204
###################################################
fls1 = BMS:::.flsresultlist('fls1')


###################################################
### code chunk number 20: bmsmanual.Rnw:207-208
###################################################
summary(fls1)


###################################################
### code chunk number 21: bmsmanual.Rnw:212-213
###################################################
plotConv(fls1)


###################################################
### code chunk number 22: bmsmanual.Rnw:217-218
###################################################
plotConv(fls1[1:100])


###################################################
### code chunk number 23: bmsmanual.Rnw:226-227
###################################################
pmp.bma(fls1)[1:5,]


###################################################
### code chunk number 24: bmsmanual.Rnw:230-231
###################################################
colSums(pmp.bma(fls1))


###################################################
### code chunk number 25: bmsmanual.Rnw:234-235
###################################################
coef(fls1)[1:5,]


###################################################
### code chunk number 26: bmsmanual.Rnw:238-239
###################################################
coef(fls1,exact=TRUE)[1:5,]


###################################################
### code chunk number 27: bmsmanual.Rnw:250-251
###################################################
fls2 = BMS:::.flsresultlist('fls2')


###################################################
### code chunk number 28: bmsmanual.Rnw:253-254 (eval = FALSE)
###################################################
## fls2= bms(datafls, burn=20000, iter=50000, g="BRIC", mprior="uniform", mcmc="rev.jump", start.value=0, user.int=F)


###################################################
### code chunk number 29: bmsmanual.Rnw:256-257
###################################################
summary(fls2)


###################################################
### code chunk number 30: bmsmanual.Rnw:261-263
###################################################
fls_combi = c(fls1,fls2)
summary(fls_combi)


###################################################
### code chunk number 31: bmsmanual.Rnw:277-278
###################################################
fls_g5 = BMS:::.flsresultlist('fls_g5')


###################################################
### code chunk number 32: bmsmanual.Rnw:280-281 (eval = FALSE)
###################################################
## fls_g5 = bms(datafls, burn=20000, iter=50000, g=5, mprior="uniform", user.int=F)


###################################################
### code chunk number 33: bmsmanual.Rnw:283-285
###################################################
coef(fls_g5)[1:5,]
summary(fls_g5)


###################################################
### code chunk number 34: bmsmanual.Rnw:300-301
###################################################
fls_ebl = BMS:::.flsresultlist('fls_ebl')


###################################################
### code chunk number 35: bmsmanual.Rnw:303-304 (eval = FALSE)
###################################################
## fls_ebl = bms(datafls, burn=20000, iter=50000, g="EBL", mprior="uniform", nmodel=1000, user.int=F)


###################################################
### code chunk number 36: bmsmanual.Rnw:306-307
###################################################
summary(fls_ebl)


###################################################
### code chunk number 37: bmsmanual.Rnw:310-311
###################################################
plot(fls_ebl)


###################################################
### code chunk number 38: bmsmanual.Rnw:320-321
###################################################
fls_hyper = BMS:::.flsresultlist('fls_hyper')


###################################################
### code chunk number 39: bmsmanual.Rnw:323-324 (eval = FALSE)
###################################################
## fls_hyper = bms(datafls, burn=20000, iter=50000, g="hyper=UIP", mprior="random", mprior.size=7, nmodel=1000, user.int=F)


###################################################
### code chunk number 40: bmsmanual.Rnw:326-327
###################################################
summary(fls_hyper)


###################################################
### code chunk number 41: bmsmanual.Rnw:331-332
###################################################
gdensity(fls_hyper)


###################################################
### code chunk number 42: bmsmanual.Rnw:336-337
###################################################
image(fls_hyper)


###################################################
### code chunk number 43: bmsmanual.Rnw:345-346
###################################################
density(fls_combi,reg="Muslim")


###################################################
### code chunk number 44: bmsmanual.Rnw:350-351
###################################################
coef(fls_combi,exact=T,condi.coef=T)["Muslim",]


###################################################
### code chunk number 45: bmsmanual.Rnw:357-358
###################################################
dmuslim=density(fls_hyper,reg="Muslim",addons="Eebl")


###################################################
### code chunk number 46: bmsmanual.Rnw:364-365
###################################################
quantile(dmuslim, c(0.025, 0.975))


###################################################
### code chunk number 47: bmsmanual.Rnw:372-375
###################################################
fcstbma= bms(datafls[1:70,], mprior="uniform", burn=20000, iter=50000, user.int=FALSE)

pdens = pred.density(fcstbma, newdata=datafls[71:72,])


###################################################
### code chunk number 48: bmsmanual.Rnw:381-382
###################################################
plot(pdens, 2)


###################################################
### code chunk number 49: bmsmanual.Rnw:388-389
###################################################
quantile(pdens, c(0.05, 0.95))


###################################################
### code chunk number 50: bmsmanual.Rnw:395-396
###################################################
pdens$dyf(datafls[71:72,1])


###################################################
### code chunk number 51: bmsmanual.Rnw:400-401
###################################################
plot(pdens, "ZM", realized.y=datafls["ZM",1])


###################################################
### code chunk number 52: bmsmanual.Rnw:408-409
###################################################
lps.bma(pdens, datafls[71:72,1])


###################################################
### code chunk number 53: bmsmanual.Rnw:498-499
###################################################
data(attitude)


###################################################
### code chunk number 54: bmsmanual.Rnw:502-503
###################################################
att_full = zlm(attitude,g="UIP")


###################################################
### code chunk number 55: bmsmanual.Rnw:506-507
###################################################
summary(att_full)


###################################################
### code chunk number 56: bmsmanual.Rnw:511-513
###################################################
att_best = as.zlm(att,model=1)
summary(att_best)


###################################################
### code chunk number 57: bmsmanual.Rnw:518-520
###################################################
att_bestlm = lm(model.frame(as.zlm(att)))
summary(att_bestlm)


###################################################
### code chunk number 58: bmsmanual.Rnw:529-530
###################################################
att_learn = bms(attitude,mprior="uniform", fixed.reg=c("complaints", "learning") )


###################################################
### code chunk number 59: bmsmanual.Rnw:539-540
###################################################
fls_culture = bms(datafls,fixed.reg=c(1,8:16,24,26:41), mprior="random", mprior.size=28, mcmc="enumeration", user.int=F)


###################################################
### code chunk number 60: bmsmanual.Rnw:544-545
###################################################
coef(fls_culture)[28:41, ]


###################################################
### code chunk number 61: bmsmanual.Rnw:549-550
###################################################
plotModelsize(fls_culture, ksubset=27:41)


