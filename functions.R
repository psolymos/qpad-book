#' Simulation based confidence nad prediction intervals for LM/GLM
predict_sim <-
function(object, newdata=NULL,
interval = c("none", "confidence", "prediction"),
type=c("asymp", "pboot", "npboot"),
level=0.95, B=99, ...) {
    interval <- match.arg(interval)
    type <- match.arg(type)
    if (is.null(newdata)) {
        x <- model.frame(object)
        X <- model.matrix(object)
    } else {
        x <- model.frame(delete.response(terms(object)), newdata)
        X <- model.matrix(attr(x, "terms"), x)
    }
    n <- nrow(x)
    fun <- switch(family(object)$family,
        "gaussian"=function(x) rnorm(length(x), x, summary(object)$sigma),
        "poisson"= function(x) rpois(length(x), x),
        "binomial"=function(x) rbinom(length(x), 1, x),
        stop("model family not recognized"))
    if (interval=="none")
        return(predict(object, newdata, ...))
    if (B < 2)
        stop("Are you kidding? B must be > 1")
    if (type == "asymp") {
        cm <- rbind(coef(object),
            MASS::mvrnorm(B, coef(object), vcov(object)))
        #fm <- apply(cm, 1, function(z) X %*% z)
    }
    if (type == "boot") {
        cm <- matrix(0, B+1, length(coef(object)))
        cm[1,] <- coef(object)
        xx <- model.frame(object)
        for (i in 2:B) {
            j <- sample.int(n, n, replace=TRUE)
            cm[i,] <- coef(update(object, data=xx[j,]))
        }
        #fm <- apply(cm, 1, function(z) X %*% z)
    }
    if (type == "npboot") {
        cm <- matrix(0, B+1, length(coef(object)))
        cm[1,] <- coef(object)
        xx <- model.frame(object)
        j <- attr(attr(xx, "terms"), "response")
        f <- fitted(object)
        for (i in 2:B) {
            xx[,j] <- fun(f)
            cm[i,] <- coef(update(object, data=xx))
        }
        #fm <- apply(cm, 1, function(z) X %*% z)
    }
    fm <- X %*% t(cm)
    fm <- family(object)$linkinv(fm)
    y <- if (interval == "prediction")
        matrix(fun(fm), n, B+1) else fm
    rownames(y) <- rownames(x)
    p <- c(0.5, (1-level) / 2, 1 - (1-level) / 2)
    stat_fun <- function(x)
        c(mean(x), sd(x), quantile(x, p))
    out <- cbind(fm[,1], t(apply(y, 1, stat_fun)))
    colnames(out) <- c("fit", "mean", "se", "median", "lwr", "upr")
    data.frame(out[,c("fit", "lwr", "upr", "mean", "median", "se")])
}
#'
#' Internal function
.r2_fun <-
function(observed, fitted, distr=c("binomial", "poisson"),
size=1, null=NULL, p=0)
{
    distr <- match.arg(distr)
    if (distr == "poisson") {
        if (is.null(null))
            null <- mean(observed)
        ll0 <- sum(dpois(observed, null, log=TRUE))
        lls <- sum(dpois(observed, observed, log=TRUE))
        llf <- sum(dpois(observed, fitted, log=TRUE))
    } else {
        if (is.null(null))
            null <- mean(observed/size)
        ll0 <- sum(dbinom(observed, size, null, log=TRUE))
        lls <- sum(dbinom(observed, size, observed/size, log=TRUE))
        llf <- sum(dbinom(observed, size, fitted, log=TRUE))
    }
    n <- length(observed)
    R2 <- 1 - (lls - llf) / (lls - ll0)
    R2adj <- 1 - (1 - R2) * ((n-1) / (n-(p+1)))
    D0 <- -2 * (ll0 - lls)
    DR <- -2 * (llf - lls)
    p_value <- 1 - pchisq(D0 - DR, p)
    #p_value <- 1 - pchisq(DR, length(observed)-(p+1))
    c(R2=R2, R2adj=R2adj, Deviance=D0 - DR, Dev0=D0, DevR=DR,
        #df=p,
        df0=length(observed)-1, dfR=length(observed)-(p+1),
        p_value=p_value)
}
#'
#' Deviance based R^2
R2dev <-
function(object, ...) {
    y <- model.response(model.frame(object), "numeric")
    f <- fitted(object)
    .r2_fun(y, fitted(object),
        distr=family(object)$family, size=1, null=NULL,
        p=length(coef(object))-1)
}
