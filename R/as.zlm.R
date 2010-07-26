as.zlm <-
function (bmao, model = 1) 
{
    thiscall = match.call()
    if (!is.bma(bmao)) 
        stop("bmao needs to be a bma object")
    bools = bmao$topmod$bool()
    if (is.character(model)) {
        model = (1:length(bools))[bools == model[[1]]]
        if (length(model) == 0) 
            stop("Provided model index was not found in bmao object topmodels")
    }
    else if ((length(model) == bmao$info$K) && (is.numeric(model) || 
        is.logical(model))) {
        model = (1:length(bools))[bools == bin2hex(model)]
        if (length(model) == 0) 
            stop("Provided model index was not found in bmao object topmodels")
    }
    else if ((length(model) == 1) && (is.numeric(model) || is.logical(model))) {
        if (model < 1 | model > length(bools)) 
            stop("Provided model index was not found in bmao object topmodels")
    }
    else stop("model needs to be an integer or other model index representation")
    inclvbls = as.logical(bmao$topmod$bool_binary()[, model, 
        drop = TRUE])
    yXdf = as.data.frame(bmao$X.data)
    zlmres = zlm(as.formula(yXdf[, c(TRUE, inclvbls)]), data = yXdf, 
        g = bmao$gprior.info)
    zlmres$call <- thiscall
    return(zlmres)
}
