bms <-
function (X.data, burn = 1000, iter = NA, nmodel = 500, mcmc = "bd", 
    g = "UIP", mprior = "random", mprior.size = NA, user.int = TRUE, 
    start.value = NA, g.stats = TRUE, logfile = FALSE, logstep = 10000, 
    force.full.ols = FALSE) 
{
    if (class(X.data)[[1]] == "formula") {
        X.data = stats::model.frame(X.data)
        if (!is.null(ncol(X.data[[2]]))) 
            X.data = cbind(X.data[[1]], X.data[[2]][, -1])
    }
    if (any(is.na(X.data))) {
        X.data = na.omit(X.data)
        if (nrow(X.data) < 3) {
            stop("Too few data observations. Please provide at least three data rows without NA entries.")
        }
        warning("Argument 'X.data' contains NAs. The corresponding rows have not been taken into account.")
    }
    N <- nrow(X.data)
    K = ncol(X.data) - 1
    maxk = N - 3
    return.g.stats = g.stats
    if (nmodel[1] <= 0 | is.na(nmodel[1])) {
        dotop = FALSE
        nmodel = 0
    }
    else {
        dotop = TRUE
    }
    if (missing(mcmc) && (K < 15)) {
        mcmc = "enum"
    }
    int = FALSE
    is.enum = FALSE
    if (length(grep("int", mcmc, ignore.case = TRUE))) {
        int = TRUE
    }
    if (length(grep("enum", mcmc, ignore.case = TRUE))) {
        is.enum = TRUE
        sampling = .iterenum
        if (K > maxk) 
            sampling = .iterenum.KgtN
    }
    else if (length(grep("bd", mcmc, ignore.case = TRUE))) {
        sampling = switch(int + 1, .fls.samp, .fls.samp.int)
    }
    else {
        sampling = switch(int + 1, .rev.jump, .rev.jump.int)
    }
    if (is.enum) {
        start.value2 = 0
        if (length(start.value) == 1) {
            start.value2 = suppressWarnings(as.integer(start.value))
            if (any(is.na(start.value2)) | start.value2[[1]] < 
                K + 5 | start.value2[[1]] < 0 | start.value2[[1]] >= 
                (2^K - 1)) {
                start.value = 0
                start.value2 = 0
            }
            else {
                start.value = .enum_fromindex(start.value2)
                start.value = c(numeric(K - length(start.value)), 
                  start.value)
            }
        }
        else {
            start.value = 0
        }
        burn = 0
        int = FALSE
        mcmc = "enum"
        is.enum = TRUE
        if (K > maxk) {
            lastindex = 2^K - 1 - sum(choose(K, (N - 2):K))
        }
        else lastindex = 2^K - 1
        if (is.na(iter)) {
            iter = lastindex - start.value2
        }
        iter = min(iter, 2^K - 1 - start.value2)
    }
    else {
        if (is.na(iter)) {
            iter = 3000
        }
    }
    if (logfile != FALSE) {
        if (is.character(logfile)) {
            sfilename = logfile
        }
        else {
            sfilename = "test.log"
        }
        if (nchar(sfilename) > 0) 
            file.create(sfilename)
        logfile = TRUE
        cat(as.character(Sys.time()), ": starting loop ... \n", 
            append = TRUE, file = sfilename)
        if (logstep != 10000) 
            fact = logstep
        else fact = max(floor((burn + iter)/100), logstep)
    }
    pmplist = .choose.mprior(mprior, mprior.size, K = K)
    mprior = pmplist$mp.mode
    y <- as.matrix(X.data[, 1])
    X <- as.matrix(X.data[, 2:ncol(X.data)])
    y.mean = mean(y)
    y <- y - matrix(y.mean, N, 1, byrow = TRUE)
    X.mean = colMeans(X)
    X <- X - matrix(X.mean, N, K, byrow = TRUE)
    XtX.big = crossprod(X)
    Xty.big = crossprod(X, y)
    yty = as.vector(crossprod(y))
    coreig = eigen(cor(X), symmetric = TRUE, only.values = TRUE)$values
    if (sum(coreig > 1e-10) < min(K, (N - 1))) {
        force.full.ols = TRUE
    }
    if (int) {
        if (length(grep("#", colnames(X.data), fixed = TRUE)) == 
            0) 
            stop("Please separate column names of interaction terms by # (e.g. A#B)")
        mPlus = .constr.intmat(X, K)
    }
    else {
        mPlus <- NA
    }
    gprior.info = .choose.gprior(g, N, K, g.stats)
    if (gprior.info$gtype == "EBL") {
        lprobcalc = .lprob.eblocal.init(N = N, K = K, yty = yty, 
            return.g = gprior.info$return.g.stats)
    }
    else if (gprior.info$gtype == "hyper") {
        lprobcalc = .lprob.hyperg.init(N = N, K = K, yty = yty, 
            f21a = gprior.info$hyper.parameter, return.gmoments = gprior.info$return.g.stats)
    }
    else {
        lprobcalc = .lprob.constg.init(g = gprior.info$g, N = N, 
            K = K, yty = yty)
    }
    start.list = .starter(K, start.value, y, N = N, XtX.big = XtX.big, 
        Xty.big = Xty.big, X = X)
    molddraw = start.list$molddraw
    start.position = start.list$start.position
    kold = sum(molddraw)
    position = (1:K)[molddraw == 1]
    collect.otherstats = FALSE
    otherstats = numeric(0)
    add.otherstats = numeric(0)
    if (gprior.info$return.g.stats & !(gprior.info$is.constant)) {
        add.otherstats = gprior.info$shrinkage.moments
        collect.otherstats = TRUE
    }
    cumsumweights = iter
    if (collect.otherstats) {
        addup = .addup.mcmc.wotherstats
    }
    else {
        addup = .addup.mcmc
    }
    if (is.enum) {
        cumsumweights = 0
        if (collect.otherstats) {
            addup = .addup.enum.wotherstats
        }
        else {
            addup = .addup.enum
        }
    }
    environment(addup) <- environment()
    ols.object = .ols.terms2(positions = (1:K)[molddraw == 1], 
        yty = yty, k = kold, N, K = K, XtX.big = XtX.big, Xty.big = Xty.big)
    lik.list = lprobcalc$lprob.all(ymy = ols.object$ymy, k = kold, 
        bhat = ols.object$bhat, diag.inverse = ols.object$diag.inverse)
    lprobold = lik.list$lprob
    b1 = lik.list$b1new
    b2 = lik.list$b2new
    pmpold = pmplist$pmp(ki = kold, mdraw = molddraw)
    null.lik = ((1 - N)/2) * log(yty)
    topmods = .top10(nmaxregressors = K, nbmodel = nmodel, bbeta = FALSE, 
        bbeta2 = FALSE, lengthfixedvec = length(add.otherstats))
    if (mcmc == "enum") {
        try(topmods$duplicates_possible(FALSE), silent = TRUE)
    }
    if (dotop) 
        topmods$addmodel(mylik = pmpold + lprobold, vec01 = molddraw, 
            fixedvec = lik.list$otherstats)
    null.count = 0
    models.visited = 0
    inccount = numeric(K)
    msize = 0
    k.vec = numeric(K)
    b1mo = numeric(K)
    ab = numeric(K)
    b2mo = numeric(K)
    bb = numeric(K)
    possign = inccount
    mnewdraw = numeric(K)
    if (force.full.ols) {
        candi.is.full.object = TRUE
    }
    else {
        candi.is.full.object = FALSE
    }
    bmo = numeric(4 * K)
    bm = bmo
    if (is.enum) {
        addup()
    }
    set.seed(as.numeric(Sys.time()))
    t1 <- Sys.time()
    nrep = burn + iter
    i = 0
    while (i < nrep) {
        i = i + 1
        if (logfile) {
            if (i%%fact == 0) {
                cat(as.character(Sys.time()), ":", i, "current draw \n", 
                  append = TRUE, file = sfilename)
            }
        }
        a = sampling(molddraw = molddraw, K = K, mPlus = mPlus, 
            maxk = maxk, oldk = kold)
        mnewdraw = a$mnewdraw
        positionnew = a$positionnew
        knew = length(positionnew)
        pmpnew = pmplist$pmp(ki = knew, mdraw = mnewdraw)
        if (!is.enum) {
            if (int) {
                if (length(c(a$dropi, a$addi)) > 2 | i < 3 | 
                  force.full.ols) {
                  candi.is.full.object = TRUE
                }
                else {
                  candi.is.full.object = FALSE
                }
            }
            if (candi.is.full.object) {
                ols.candidate = .ols.terms2(positions = positionnew, 
                  yty = yty, k = knew, N, K = K, XtX.big = XtX.big, 
                  Xty.big = Xty.big)
                ymy.candi = ols.candidate$ymy
            }
            else {
                ymy.candi = ols.object$child.ymy(a$addi, a$dropi, 
                  k = knew)
            }
            lprobnew = lprobcalc$just.loglik(ymy.candi, knew)
            accept.candi = as.logical(log(.Internal(runif(1, 
                0, 1))) < lprobnew - lprobold + pmpnew - pmpold)
        }
        else {
            accept.candi = TRUE
            candi.is.full.object = FALSE
        }
        if (accept.candi) {
            if (!candi.is.full.object) {
                ols.res = ols.object$mutate(addix = a$addi, dropix = a$dropi, 
                  newpos = positionnew, newk = knew)
            }
            else {
                ols.object = ols.candidate
                ols.res = ols.candidate$full.results()
            }
            lik.list = lprobcalc$lprob.all(ols.res$ymy, knew, 
                ols.res$bhat, ols.res$diag.inverse)
            lprobold = lik.list$lprob
            position = positionnew
            pmpold = pmpnew
            molddraw = mnewdraw
            kold = knew
            models.visited = models.visited + 1
        }
        if (i > burn) {
            b1 = lik.list$b1new
            b2 = lik.list$b2new
            addup()
            if (dotop) 
                topmods$addmodel(mylik = pmpold + lprobold, vec01 = molddraw, 
                  fixedvec = otherstats)
        }
    }
    if (dotop) 
        topmods = .topmod.as.bbetaT(topmods, gprior.info, X.data)
    timed <- difftime(Sys.time(), t1)
    if (is.enum) {
        iter = iter + 1
        models.visited = models.visited + 1
    }
    bmo = matrix(bmo, 4, byrow = TRUE)
    b1mo = bmo[1, ]
    b2mo = bmo[2, ]
    k.vec = bmo[3, ]
    possign = bmo[4, ]
    rm(bmo)
    post.inf = .post.calc(gprior.info, add.otherstats, k.vec, 
        null.count, X.data, topmods, b1mo, b2mo, iter, burn, 
        inccount, models.visited, K, N, msize, timed, cumsumweights, 
        mcmc, possign)
    result = list(info = post.inf$info, arguments = .construct.arglist(bms), 
        topmod = topmods, start.pos = sort(start.position), gprior.info = post.inf$gprior.info, 
        mprior.info = pmplist, X.data = X.data, reg.names = post.inf$reg.names, 
        bms.call = match.call(bms, sys.call(0)))
    class(result) = c("bma")
    if (user.int) {
        print(result)
        print(timed)
        plot.bma(result)
    }
    return(invisible(result))
}
