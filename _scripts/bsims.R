#' # bSims: bird simulations
#'
#' ## Landscape initialization
#'
library(intrval)
library(MASS)
library(ADPclust)
library(mefa4)
#library(spatstat)
#'
#'


bsims_init <- function(
  extent=10,
  road=0,
  edge=0,
  offset=0
) {
  if (extent <= 0)
    stop("extent must be positive")
  road <- abs(road)
  edge <- abs(edge)
  a <- extent / 2
  box <- cbind(
    x=rep(a, 4) * c(-1, -1, 1, 1),
    y=rep(a, 4) * c(-1, 1, 1, -1))
  strata <- c(
    -a,
    offset-road-edge,
    offset-road,
    offset+road,
    offset+road+edge,
    a)
  strata <- pmin(a, pmax(-a, strata))
  names(strata) <- c("+h", "he", "er", "re", "eh", "h+")
  x <- list(
    extent=extent,
    road=road,
    edge=edge,
    offset=offset,
    box=box,
    strata=strata)
  class(x) <- c("bsim", "bsims_landscape")
  x
}

print.bsims_landscape <- function(x, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  her <- paste0(
    ifelse(A[1] > 0, "H", ""),
    ifelse(A[2] > 0, "E", ""),
    ifelse(A[3] > 0, "R", ""), collapse="")
  cat("bSims landscape\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her, "\n", sep="")
  invisible(x)
}

plot.bsims_landscape <-
function(x,
  col=c("darkolivegreen1", "burlywood1", "lightgrey"),
  xlim=NULL, ylim=NULL, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  plot(0, type="n",
    xlim=if (is.null(xlim)) range(x$box[,"x"]) else xlim,
    ylim=if (is.null(ylim)) range(x$box[,"y"]) else ylim,
    xlab="", ylab="", axes=FALSE, asp=1, ...)
  if (A[1] > 0)
    polygon(x$box, col=col[1], border=NA)
  if (A[2] > 0)
    polygon(x$strata[c("he", "he", "eh", "eh")],
      x$strata[c("+h", "h+", "h+", "+h")],
      col=col[2], border=NA)
  if (A[3] > 0)
    polygon(x$strata[c("er", "er", "re", "re")],
      x$strata[c("+h", "h+", "h+", "+h")],
      col=col[3], border=NA)
  invisible(x)
}


#' ## Population
#' Can be extended to modify uniform/uniform into clustered, systematic, etc
#' to represent deviation from random
acceptreject <- function(
  n, # number of pts
  f=NULL, # function returning probability given distance
  x0=0, x1=1, # x range
  y0=0, y1=1, # y range
  m=0, # margin width for edge effect
  maxit=100, # max iterations to break
  fail=FALSE
) {
  if (n <= 0)
    return(cbind(
        x=numeric(0),
        y=numeric(0)))
  if (is.null(f))
      return(cbind(
        x=runif(n, x0, x1),
        y=runif(n, y0, y1)))
  g <- function(m)
    c(runif(1, x0-m, x1+m), runif(1, y0-1, y1+m))
  xy <- rbind(g(m=0))
  colnames(xy) <- c("x", "y")
  i <- 1L # iterations
  j <- 1L # pts inside bbox
  while (j < n) {
    if (i > maxit*n)
      break
    pr <- g(m=m)
    d <- min(sqrt((xy[,1]-pr[1])^2 + (xy[,2]-pr[2])^2))
    p <- f(d)
    r <- runif(1)
    if (r <= p)
      xy <- rbind(xy, unname(pr))
    i <- i + 1L
    j <- sum(xy[,1] %[]% c(x0, x1) & xy[,2] %[]% c(y0, y1))
  }
  xy <- xy[xy[,1] %[]% c(x0, x1) & xy[,2] %[]% c(y0, y1),,drop=FALSE]
  if (nrow(xy) < n) {
    msg <- sprintf("number of points generated: %.0f", nrow(xy))
    if (fail)
      stop(msg)
    warning(msg)
    xy <- rbind(xy, matrix(NA, n-nrow(xy), 2))
  }
  xy
}

bsims_populate <- function(
  x, # landscape object
  density=1, # D, recycled 3x for HER
  abund_fun=NULL, # N ~ Pois(lambda), lambda=DA
  xy_fun=NULL, # NULL ~ CSR complete spatial randomness
  margin=0, # margin to pass to xy_fun for edge effect, units as in extent
  maxit=100, # x N times to try
  fail=FALSE,
  ...)
{
  if (!inherits(x, "bsims_landscape"))
    stop("x must be a bsims_landscape object")
  A <- diff(x$strata) * diff(range(x$strata))
  D <- rep(density, 3)[c(1,2,3,2,1)]
  lambda <- A * D
  if (is.null(abund_fun))
    abund_fun <- function(lambda, ...) rpois(1, lambda)
  N <- sapply(lambda, abund_fun, ...)
  names(A) <- names(D) <- names(N) <- names(lambda) <-
    c("+H", "+E", "R", "E+", "H+")
  d <- NULL
  for (i in 1:5) {
    z <-data.frame(
      acceptreject(n=N[i], f=xy_fun,
        x0=x$strata[i], x1=x$strata[i+1],
        y0=x$strata[1], y1=x$strata[6],
        m=margin,
        maxit=maxit,
        fail=fail),
      s=rep(i, N[i]))
    # add here spatial non-randomness
    d <- rbind(d, z)
  }
  d$i <- seq_len(nrow(d))
  d$s <- factor(c("H", "E", "R", "E", "H")[d$s], c("H", "E", "R"))
  d <- d[,c("i", "s", "x", "y")]
  x$nests <- d
  x$abundance <- N
  x$lambda <- lambda
  x$area <- A
  x$density <- D
  class(x) <- c("bsim", "bsims_population")
  x
}

print.bsims_population <- function(x, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  her <- paste0(
    ifelse(A[1] > 0, "H", ""),
    ifelse(A[2] > 0, "E", ""),
    ifelse(A[3] > 0, "R", ""), collapse="")
  cat("bSims population\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her,
    "\n  total abunance: ", sum(x$abundance), "\n", sep="")
  invisible(x)
}

points.bsims_population <-
function(x, ...) {
  points(x$nests[,c("x", "y")], ...)
  invisible(x)
}

plot.bsims_population <-
function(x, pch_nest=3,
col_nest="darkgreen", cex_nest=1,
...) {
  op <- par(xpd = TRUE)
  on.exit(par(op))
  xx <- x
  class(xx) <- c("bsim", "bsims_landscape")
  plot(xx, ...)
  if (!is.na(pch_nest))
    points(x, pch=pch_nest, col=col_nest, cex=cex_nest, ...)
  invisible(x)
}

timetoevent <- function(rate, duration) {
  te <- rexp(n=duration*rate, rate=rate)
  while(sum(te) < duration) {
    te <- c(te, rexp(n=duration*rate, rate=rate))
  }
  cte <- cumsum(te)
  te[cte < duration]
}
# avoid can be used to limit the movement in x direction
# it gives an interval to avoid
# note: nest location must be accounted for
events <-
function(vocal_rate=1, move_rate=1,
duration=10, movement=0, avoid=c(0,0)) {
  ev <- cumsum(timetoevent(vocal_rate, duration))
  em <- cumsum(timetoevent(move_rate, duration))
  iv <- rep(1, length(ev))
  im <- rep(0, length(em))
  dv <- matrix(NA, length(ev), 2)
  dm <- MASS::mvrnorm(length(em), c(0, 0), diag(movement^2, 2, 2))
  dm <- dm[dm[,1] %][% avoid,,drop=FALSE]
  while (nrow(dm) < length(em)) {
    dm <- rbind(dm, MASS::mvrnorm(length(em), c(0, 0), diag(movement^2, 2, 2)))
    dm <- dm[dm[,1] %][% avoid,,drop=FALSE]
  }
  dm <- dm[seq_along(em),,drop=FALSE]
  ## don't store coordinates that are same as previous for movement
  if (nrow(dm) > 1L) {
    keep <- logical(nrow(dm))
    keep[1L] <- TRUE
    for (i in 2:length(keep)) {
      keep[i] <- !all(dm[i,] == dm[i-1L,])
    }
  }
  h <- cbind(rbind(dv, dm), c(ev, em), c(iv, im))
  colnames(h) <- c("x", "y", "t", "v")
  o <- order(h[,"t"])
  h <- as.data.frame(h[o,,drop=FALSE])
  for (i in which(is.na(h[,"x"]))) {
    if (i == 1L) {
      h$x[i] <- 0
      h$y[i] <- 0
    } else {
      h$x[i] <- h$x[i-1L]
      h$y[i] <- h$y[i-1L]
    }
  }
  h
}

## rate can be
## - single number
## - vector matching length of mixture
## - vector of length 3 with mixture=1: ~ HER
## - matrix of 3 x mixture: ~ HER + groups
bsims_animate <- function(
  x, # population object
  vocal_rate=1, # phi /min
  move_rate=1, #movement
  duration=10,
  movement=0, # SD for 2D kernel
  mixture=1, # finite mixture group proportions
  avoid=c("none", "R", "ER"),
  ...) {
  if (!inherits(x, "bsims_population"))
    stop("x must be a bsims_population object")
  avoid <- match.arg(avoid)
  if (avoid == "ER" && sum(x$density[2:4]) > 0)
    stop(">0 density ER strata cannot be avoided")
  if (avoid == "R" && sum(x$density[3]) > 0)
    stop(">0 density R stratum cannot be avoided")
  if (any(mixture < 0))
    stop("mixture must not be negative")
  K <- length(mixture)
  G <- paste0("G", 1:K)
  P <- structure(mixture / sum(mixture), names=G)
  ## vocal rate processing
  if (length(vocal_rate) == 1L) {
    vr <- matrix(vocal_rate, 3, K)
  } else {
    if (is.null(dim(vocal_rate))) {
      if (K > 1L) {
        if (length(vocal_rate) != K)
          stop("vocal_rate length must equal mixture length")
        vr <- matrix(vocal_rate, 3, K, byrow=TRUE)
      } else {
        if (length(vocal_rate) != 3L)
          stop("vocal_rate length must equal 3 when length(mixture)=1")
        vr <- matrix(vocal_rate, 3, K)
      }
    } else {
      if (dim(vocal_rate) != c(3L, K))
        stop("vocal_rate dimension must be 3 x length(mixture)")
      vr <- vocal_rate
    }
  }
  ## movement rate processing
  if (length(move_rate) == 1L) {
    mr <- matrix(move_rate, 3, K)
  } else {
    if (is.null(dim(move_rate))) {
      if (K > 1L) {
        if (length(move_rate) != K)
          stop("move_rate length must equal mixture length")
        mr <- matrix(move_rate, 3, K, byrow=TRUE)
      } else {
        if (length(move_rate) != 3L)
          stop("move_rate length must equal 3 when length(mixture)=1")
        mr <- matrix(move_rate, 3, K)
      }
    } else {
      if (dim(move_rate) != c(3L, K))
        stop("move_rate dimension must be 3 x length(mixture)")
      mr <- move_rate
    }
  }
  dimnames(vr) <- dimnames(mr) <- list(c("H", "E", "R"), G)
  N <- sum(x$abundance)
  g <- sample(G, N, replace=TRUE, prob=P)
  x$nests$g <- factor(g, G)
  s <- as.character(x$nests$s)
  Events <- list()
  for (i in seq_len(N)) {
    a <- switch(avoid,
      "none" = c(0,0),
      "R" = x$strata[c("er", "re")]-x$nests$x[i],
      "ER" = x$strata[c("he", "eh")]-x$nests$x[i])
    Events[[i]] <- events(
      vocal_rate=vr[s[i], g[i]],
      move_rate=mr[s[i], g[i]],
      duration=duration, movement=movement, avoid=a)
  }
  x$vocal_rate <- vr
  x$move_rate <- mr
  x$duration <- duration
  x$movement <- movement
  x$mixture <- P
  x$events <- Events
  class(x) <- c("bsim", "bsims_events")
  x
}
print.bsims_events <- function(x, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  her <- paste0(
    ifelse(A[1] > 0, "H", ""),
    ifelse(A[2] > 0, "E", ""),
    ifelse(A[3] > 0, "R", ""), collapse="")
  cat("bSims events\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her,
    "\n  total abunance: ", sum(x$abundance),
    "\n  ", ifelse(length(x$mixture) > 0, "mixture with ", ""),
    "total duration: ", x$duration, "\n", sep="")
  invisible(x)
}

get_events <- function(x, vocal_only=TRUE) {
  if (sum(x$abundance) == 0)
    return(data.frame(
      x=numeric(0),
      y=numeric(0),
      t=numeric(0),
      v=numeric(0),
      i=numeric(0)
    ))
  z <- lapply(1:length(x$events), function(i) {
    zz <- x$events[[i]]
    zz$i <- i
    zz
  })
  z <- do.call(rbind, z)
  z <- z[order(z$t),]
  if (vocal_only)
    z <- z[z$v > 0,,drop=FALSE]
  rownames(z) <- NULL
  z$x <- x$nests$x[z$i] + z$x
  z$y <- x$nests$y[z$i] + z$y
  z
}
points.bsims_events <-
function(x, vocal_only=TRUE, ...) {
  points(get_events(x, vocal_only), ...)
  invisible(x)
}
lines.bsims_events <-
function(x, ...) {
  N <- length(x$events)
  for (i in seq_len(N)) {
    xy <- cbind(
      x=x$nests$x[i] + x$events[[i]]$x,
      y=x$nests$y[i] + x$events[[i]]$y)
      lines(xy, ...)
  }
  invisible(x)
}
plot.bsims_events <-
function(x,
pch_nest=3, col_nest="darkgreen", cex_nest=1,
pch_vocal=21, col_vocal="blue", cex_vocal=0.5,
lty_move=1, col_move="orange", lwd_move=1,
...) {
  op <- par(xpd = TRUE)
  on.exit(par(op))
  xx <- x
  class(xx) <- c("bsim", "bsims_population")
  plot(xx, pch_nest=pch_nest, col_nest=col_nest, cex_nest=cex_nest, ...)
  if (!is.na(lty_move))
    lines(x, col=col_move, lty=lty_move, lwd=lwd_move, ...)
  if (!is.na(pch_vocal))
    points(x, vocal_only=TRUE,
      col=col_vocal, pch=pch_vocal, cex=cex_vocal, ...)
  invisible(x)
}

## detect
## HER case need to deal with sound attenuation

## d is distance
## tau is single param for dist_fun, need to be in order
## b: breaks for stratum boundaries, 1 less in length than tau
##    must be positive
dist_fun2 <- function(d, tau, dist_fun, b=numeric(0)) {
  b <- sort(b)
  h <- 1
  for (i in seq_len(length(b))) {
    h <- c(h,
      h[length(h)] * dist_fun(b[i], tau[i]) /
        dist_fun(b[i], tau[i+1L]))
  }
  j <- cut(d, c(0, b, Inf), labels=FALSE, include.lowest=TRUE)
  dist_fun(d, tau[j]) * h[j]
}
if (FALSE) {
tau <- c(1, 2, 3, 2, 1)
tauseq <- 3:1 # c(1,2,3) # HER=123 sequence for tau values
d <- 0:400/100
plot(d, dist_fun2(d, tau[1], dist_fun), type="l")
lines(d, dist_fun2(d, tau[2], dist_fun))
lines(d, dist_fun2(d, tau[3], dist_fun))
b <- c(0.5, 1, 1.5, 2) # points at the HER boundaries
abline(v=b)
lines(d, dist_fun2(d, tau, dist_fun, b), col=2, lwd=3)
}
bsims_detect <- function(
  x,
  xy=c(0,0), # observer location
  tau=1, # can vector when HER attenuation used, compatible w/ dist_fun
  dist_fun=NULL, # takes args d and tau (single parameter)
  repel=0, # radius within which vocalizations are invalidated
  ...)
{
  if (!inherits(x, "bsims_events"))
    stop("x must be a bsims_events object")
  xy <- as.numeric(xy[1:2])
  if (any(xy %)(% range(x$strata)))
    stop("observer xy must be within extent")
  if (is.null(dist_fun))
    dist_fun <- function(d, tau) exp(-d^2/tau^2)
  N <- sum(x$abundance)
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  ## need for *att*enuation
  att <- length(tau) > 1 && A["H"] < sum(A)
  if (att && length(tau) != 3L)
    stop("tau length must be 3 for HER attenuation")
  if (att) {
    st <- x$strata
    st[1L] <- -Inf
    st[6L] <- Inf
    #ST <- c("H", "E", "R", "E", "H")
    sobs <- cut(xy[1L], st, labels=FALSE)
  }
  for (i in seq_len(N)) {
    z <- x$events[[i]]
    xx <- x$nests$x[i] + z$x - xy[1L]
    yy <- x$nests$y[i] + z$y - xy[2L]
    z$d <- sqrt((xx)^2 + (yy)^2)
    ## repel inds within repel distance (0 just for temp object z)
    z$v[z$d < repel] <- 0
    ## NA is placeholder for these vocalizations in return object
    x$events[[i]]$v[z$d < repel] <- NA
    keep <- z$v > 0
    z <- z[keep,,drop=FALSE]
    xx <- xx[keep]
    yy <- yy[keep]
    ## angle in degrees counter clockwise from x axis
#    a <- 180 * atan2(yy, xx) / pi
#    a[a < 0] <- 360+a[a < 0]
#    z$a <- a
    ## this is where HER attenuation somes in
    if (att) {
      theta <- atan2(yy, xx) # angle in rad
      sbrd <- cut(xx, st, labels=FALSE)
      for (j in seq_len(nrow(z))) {
        ## order tau as HEREH
        TAU <- tau[c(1,2,3,2,1)][sobs:sbrd[j]]
        ## calculate distance breaks from x and theta
        if (length(TAU) == 1) {
          b <- numeric(0)
        } else {
          ## this gives breaks along x axis
          if (sobs < sbrd[j]) { # bird right of observer
            stj <- st[(sobs+1):sbrd[j]]
          } else { # bird left of observer
            stj <- st[sobs:(sbrd[j]+1)]
          }
          ## breaks as radial distance: r=x/cos(theta)
          b <- (stj - xy[1L]) / cos(theta[j])
        }
        ## calculate q
        q <- dist_fun2(z$d[j], TAU, dist_fun, b)
      }
    } else {
      q <- dist_fun(z$d, tau)
    }
    u <- runif(length(z$d))
    z$det <- ifelse(u <= q, 1, 0) # detected
    z <- z[z$det > 0,,drop=FALSE]
    ## error is shown where detected, NA when not detected
    x$events[[i]]$d <- z$d[match(rownames(x$events[[i]]), rownames(z))]
  }
  x$xy <- xy
  x$tau <- tau
  x$repel <- repel
  class(x) <- c("bsim", "bsims_detections")
  x
}
print.bsims_detections <- function(x, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  her <- paste0(
    ifelse(A[1] > 0, "H", ""),
    ifelse(A[2] > 0, "E", ""),
    ifelse(A[3] > 0, "R", ""), collapse="")
  ndet <- if (sum(x$abundance) == 0)
    0 else sum(sapply(x$events, function(z) any(!is.na(z$d))))
  cat("bSims detections\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her,
    "\n  total abunance: ", sum(x$abundance),
    "\n  ", ifelse(length(x$mixture) > 0, "mixture with ", ""),
    "total duration: ", x$duration, "\n  detected: ", ndet, "\n", sep="")
  invisible(x)
}
## this adds next xy to movement xy
get_detections <- function(x, first_only=TRUE) {
  if (sum(x$abundance) == 0)
    return(data.frame(
      x=numeric(0),
      y=numeric(0),
      t=numeric(0),
      v=numeric(0),
      d=numeric(0),
      i=numeric(0)
    ))
  z <- lapply(1:length(x$events), function(i) {
    zz <- x$events[[i]]
    zz$i <- i
    zz <- zz[!is.na(zz$d),,drop=FALSE]
    zz
  })
  z <- do.call(rbind, z)
  z <- z[order(z$t),]
  if (first_only)
    z <- z[!duplicated(z$i),,drop=FALSE]
  rownames(z) <- NULL
  z$x <- x$nests$x[z$i] + z$x
  z$y <- x$nests$y[z$i] + z$y
  attr(z, "observer") <- x$xy
  z
}
points.bsims_detections <-
function(x, first_only=TRUE, ...) {
  points(get_detections(x, first_only)[,c("x", "y")], ...)
  invisible(x)
}
lines.bsims_detections <-
function(x, first_only=TRUE, ...) {
  xy <- get_detections(x, first_only)[,c("x", "y")]
  segments(
    x0=rep(x$xy[1L], nrow(xy)),
    y0=rep(x$xy[2L], nrow(xy)),
    x1=xy[,1L], y1=xy[,2L], ...)
  invisible(x)
}

plot.bsims_detections <-
function(x, first_only=TRUE,
pch_nest=3, col_nest="darkgreen", cex_nest=1,
pch_vocal=21, col_vocal="blue", cex_vocal=0.5,
lty_move=1, col_move="orange", lwd_move=1,
lty_det=1, col_det="black", lwd_det=1,
...) {
  op <- par(xpd = TRUE)
  on.exit(par(op))
  xx <- x
  class(xx) <- c("bsim", "bsims_events")
  plot(xx,
    pch_nest=pch_nest, col_nest=col_nest, cex_nest=cex_nest,
    pch_vocal=pch_vocal, col_vocal=col_vocal, cex_vocal=cex_vocal,
    lty_move=lty_move, col_move=col_move, lwd_move=lwd_move, ...)
  if (!is.na(lty_det))
    lines(x, first_only, col=col_det, lty=lty_det, ...)
  invisible(x)
}


## lognormal parametrized as mean (ybar) and SDlog
rlnorm2 <- function(n, mean = exp(0.5), sdlog = 1) {
  rlnorm(n, log(mean) - sdlog^2/2, sdlog)
}
#summary(rlnorm2(10^6, 1.3, 0.5))

bsims_transcribe <- function(
  x,
  rint=Inf,
  tint=NULL,
  first_only=TRUE,
  error=0,
  ...) {
  if (!inherits(x, "bsims_detections"))
    stop("x must be a bsims_detections object")
  tint <- if (is.null(tint))
    x$duration else sort(tint)
  if (any(tint <= 0))
    stop("tint must be > 0")
  if (any(tint > x$duration))
    stop("tint must <= duration")
  rint <- sort(rint)
  if (any(rint <= 0))
    stop("rint must be > 0")
  det <- get_detections(x, first_only)
  det <- det[det$d <= max(rint),,drop=FALSE]
  if (error < 0)
    stop("error must be >= 0")
  derr <- if (error > 0)
    rlnorm2(nrow(det), det$d, error) else det$d
  det$error <- derr - det$d
  rLAB <- paste0(c(0, round(100*rint[-length(rint)])),
    ifelse(is.finite(rint), paste0("-", round(100*rint)), "+"), "m")
  tLAB <- paste0(c(0, round(tint[-length(tint)], 2)), "-", tint, "min")
  det$rint <- factor(rLAB[cut(derr, c(0, rint), labels=FALSE,
    include.lowest=TRUE)], rLAB)
  det$tint <- factor(tLAB[cut(det$t, c(0, tint), labels=FALSE,
    include.lowest=TRUE)], tLAB)
  xt <- as.matrix(Xtab(~ rint + tint, det))
  x$detections <- det
  x$counts <- xt
  x$rint <- rint
  x$tint <- tint
  x$first_only <- first_only
  x$error <- error
  class(x) <- c("bsim", "bsims_transcript")
  x
}

## spatial patterns
## stupid
f <- function(d) ifelse(d > 0, 0, 0)
acceptreject(2, f)
try(acceptreject(2, f, fail=TRUE))
acceptreject(2, NULL)


## random
f <- function(d) ifelse(d > 0, 1, 1)
plot(seq(0,1,0.01), f(seq(0,1,0.01)), type="l")
plot(acceptreject(100, f))
nrow(acceptreject(10, f))

## systematic
f <- function(d) 1-exp(-d^2/0.1^2)
plot(seq(0,1,0.01), f(seq(0,1,0.01)), type="l")
plot(acceptreject(100, f, m=1))

## bimodal/clustered
f <- function(d) pmax(ifelse(d < 0.1, 1, 0), 0.5*(1-exp(-d^2/0.5^2)))
plot(seq(0,1,0.01), f(seq(0,1,0.01)), type="l")
plot(acceptreject(100, f, m=1))



(bsims_init())
(l <- bsims_init(3, 0.05, 0.05))
plot(l)
x <- bsims_populate(l, c(2,1,0))
plot(x)

plot(bsims_populate(l, 10))

## systematic
f <- function(d) 1-exp(-d^2/0.3^2)
plot(seq(0,1,0.01), f(seq(0,1,0.01)), type="l")
plot(x <- bsims_populate(l, 10, xy_fun=f, margin=1))

## clustered -- need to watch maxit
f <- function(d) pmax(ifelse(d < 0.2, 1, 0), 0.5*(1-exp(-d^2/1^2)))
plot(seq(0,1,0.01), f(seq(0,1,0.01)), type="l")
plot(x <- bsims_populate(l, 10, xy_fun=f, margin=1))


dim(x$nests)
sum(x$abundance)


library(detect)
l <- bsims_init(10)
p <- bsims_populate(l, 10)
a <- bsims_animate(p, phi=1) # this fails with Error in dm[, 1] : incorrect number of dimensions
o <- bsims_detect(a, tau=1)
x <- bsims_transcribe(o, rint=c(0.5, 1, Inf), tint=c(3,5,10))
plot(o, pch_vocal=NA)
x$counts
Y1 <- matrix(colSums(x$counts), nrow=1)
D1 <- matrix(x$tint, nrow=1)
Y2 <- matrix(rowSums(x$counts), nrow=1)
D2 <- matrix(x$rint, nrow=1)
exp(cmulti.fit(Y1, D1, type="rem")$coef)
exp(cmulti.fit(Y2, D2, type="dis")$coef)

## avoid R stratum
l <- bsims_init(10, 0.5, 0.5)
plot(l)
p <- bsims_populate(l, c(1,1,0))
plot(p)
a <- bsims_animate(p, movement=0.2, avoid="R")
plot(a)
o <- bsims_detect(a, c(0.5,0), repel=1)
plot(o)
## zoom in
plot(o, xlim=c(-3,3), ylim=c(-3,3))

## repel birds
l <- bsims_init(4, 0.5, 0.5)
p <- bsims_populate(l, 1)
a <- bsims_animate(p, movement=0)
o <- bsims_detect(a, c(0,0), repel=1)
plot(o)

## roadside EDR
l <- bsims_init(10, 0.5, 0.5)
p <- bsims_populate(l, 3)
a <- bsims_animate(p, movement=0)
o <- bsims_detect(a, c(0,0), tau=1:3)
plot(o, pch.point=NA)

library(magrittr)

p <- bsims_init(3) %>%
  bsims_populate(c(2,1,0)) %>%
  bsims_animate(movement=0.1)
plot(p)

rr <- 1
tt <- timetoevent(rr, 10)
op <- par(mfrow=c(1,2))
plot(ecdf(tt))
curve(1-exp(-rr*x), add=TRUE, col=2)

plot(stepfun(sort(tt), 0:length(tt)/length(tt)))
curve(1-exp(-rr*x), add=TRUE, col=2)
par(op)

## get coords
xy <- do.call(rbind, lapply(1:length(x$events), function(i) {
  cbind(x$events[[i]]$x+x$nests$x[i],x$events[[i]]$y+x$nests$y[i])
}))
## get individual id
i <- do.call(c, lapply(1:length(x$events), function(i)
  rep(i, nrow(x$events[[i]]))))
## fit ADP clustering
library(ADPclust)
ad <- adpclust(xy, dmethod = "euclidean")
## number of clusters found
ad$nclust
## classification
tab <- table(inds=i, clust=ad$clusters)




