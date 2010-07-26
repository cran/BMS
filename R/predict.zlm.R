predict.zlm <-
function (object, newdata = NULL, ...) 
{
    if (!any(class(object) == "zlm")) {
        stop("you need to provide a zlm object")
        return()
    }
    betas = object$coefficients[-1, drop = FALSE]
    alpha = object$coefficients[[1]]
    if (is.null(newdata)) {
        newX <- as.matrix(object$model[, -1, drop = FALSE])
    }
    else {
        newX = as.matrix(newdata)
        if (!is.numeric(newX)) 
            stop("newdata must be numeric!")
        if (is.vector(newdata)) 
            newX = matrix(newdata, 1)
        if (ncol(newX) != length(betas)) {
            if (ncol(newX) == length(betas) + 1) {
                newX = newX[, -1, drop = FALSE]
            }
            else {
                stop("newdata must be a matrix or data.frame with", 
                  length(betas), "columns.")
            }
        }
    }
    return(as.vector(newX %*% betas) + alpha)
}
