plotComp <-
function (..., varNr = NULL, comp = "PIP", exact = FALSE, cex.axis = 0.45, 
    main = NULL, type = "p", lty = 1:5, lwd = 1.5, pch = NULL, 
    col = NULL, cex = NULL, bg = NA, xlab = "", ylab = NULL) 
{
    bmaList = list(...)
    bmaNr = length(bmaList)
    if (bmaNr < 2) {
        stop("Submit at least two bma objects to compare results")
    }
    oldpar = par()
    if (is.null(col)) {
        col = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
            "#E6AB02", "#A6761D", "#666666")
    }
    xMat = lapply(bmaList, function(x) rownames(estimates.bma(x, 
        exact = exact)))
    xNames = xMat[[1]]
    ind = as.numeric(unlist(lapply(xMat, function(x) length(x))))
    if (length(unique(ind) > 1)) {
        smallestSet = which.min(ind)
        indMat = array(0:0, dim = c(length(xMat[[smallestSet]]), 
            length(ind)))
        for (i in 1:length(ind)) {
            indMat[, i] = as.numeric(xMat[[smallestSet]] %in% 
                xMat[[i]])
        }
        xNamesInd = which(rowSums(indMat) == bmaNr)
        xNames = xMat[[smallestSet]][xNamesInd]
    }
    compNames = c(colnames(estimates.bma(bmaList[[1]])), "Std Mean", 
        "Std Coef")
    if (is.null(xNames)) {
        stop("the bma objects have to have (the same) rownames attached to them")
    }
    if (!(comp %in% compNames)) {
        stop("Please specify comp as one of PIP, Post Mean, Post SD, Std Mean, or Std Coef")
    }
    if (comp == "Std Mean") {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            std.coefs = TRUE, exact = exact)[xNames, "Post Mean"])
        comp = "Standardized Coefficients"
    }
    else if (comp == "Std SD") {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            std.coefs = TRUE, exact = exact)[xNames, "Post SD"])
        comp = "Standardized SD"
    }
    else {
        compMatrix = sapply(bmaList, function(x) estimates.bma(x, 
            exact = exact)[xNames, comp])
    }
    bmaNames = names(bmaList)
    if (!is.null(bmaNames) & length(bmaNames) == ncol(compMatrix)) {
        colnames(compMatrix) <- bmaNames
    }
    else {
        colnames(compMatrix) = paste("Model", 1:bmaNr)
    }
    if (!is.null(varNr)) {
        compMatrix = compMatrix[varNr, , drop = FALSE]
    }
    oldmar = par()$mar
    par(mar = c(6, 4, 4, 2))
    if (is.null(ylab)) 
        ylab = paste(comp)
    if (is.null(pch)) 
        pch = 1:bmaNr
    matplot(compMatrix, main = main, type = type, col = col, 
        cex = cex, bg = bg, ylab = ylab, xlab = xlab, xaxt = "n", 
        pch = pch, lty = lty, lwd = lwd)
    legend("topright", colnames(compMatrix), pch = 1:bmaNr, col = col, 
        bty = "n")
    grid()
    axis(1, las = 2, at = 1:nrow(compMatrix), label = rownames(compMatrix), 
        cex.axis = cex.axis)
    layout(matrix(1))
    par(mar = oldmar)
}
