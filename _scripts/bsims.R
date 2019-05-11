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
  if (rate == Inf)
    return(rexp(1, rate))
  if (rate == 0)
    rate <- .Machine$double.eps
  te <- rexp(n=ceiling(duration*rate), rate=rate)
  while(sum(te) < duration) {
    te <- c(te, rexp(n=ceiling(duration*rate), rate=rate))
  }
  cte <- cumsum(te) < duration
#  if (sum(cte) < 1)
#    te[1] else te[cte]
  te[cte]
}
# avoid can be used to limit the movement in x direction
# it gives an interval to avoid
# note: nest location must be accounted for
# rmvn is a shim of MASS:mvrnorm to deal with n=1 case to return matrix
rmvn <- function(n=1L, mu, Sigma, ...) {
  if (n < 2L) {
    if (n == 0L) {
      out <- matrix(numeric(0), nrow=0L, ncol=length(mu))
    } else {
      out <- matrix(MASS::mvrnorm(n, mu, Sigma, ...), nrow=1L)
    }
    colnames(out) <- names(mu)
  } else {
    out <- MASS::mvrnorm(n, mu, Sigma, ...)
  }
  out
}
#rmvn(0, c(a=0, b=0), diag(1, 2, 2))
#rmvn(1, c(a=0, b=0), diag(1, 2, 2))
#rmvn(2, c(a=0, b=0), diag(1, 2, 2))

## if there is no movement, there is no point in making a move
## but we might want to deal with visual cues independent of movement
## so let's just leave it there
## move_rate will lead to repeated events, so that heard/seen is meaningful
## set move_rate to Inf
events <-
function(vocal_rate=1, move_rate=1,
duration=10, movement=0, avoid=c(0,0)) {
  ev <- cumsum(timetoevent(vocal_rate, duration))
  em <- cumsum(timetoevent(move_rate, duration))
  iv <- rep(1, length(ev))
  im <- rep(0, length(em))
  dv <- matrix(NA, length(ev), 2)
  dm <- rmvn(length(em), c(0, 0), diag(movement^2, 2, 2))

  dm <- dm[dm[,1] %][% avoid,,drop=FALSE]
  while (nrow(dm) < length(em)) {
    dm <- rbind(dm, rmvn(length(em), c(0, 0), diag(movement^2, 2, 2)))
    dm <- dm[dm[,1] %][% avoid,,drop=FALSE]
  }
  dm <- dm[seq_along(em),,drop=FALSE]
  h <- cbind(rbind(dv, dm), c(ev, em), c(iv, im))
  colnames(h) <- c("x", "y", "t", "v")
  o <- order(h[,"t"])
  h <- as.data.frame(h[o,,drop=FALSE])
  ## take previous known location for vocalizations
  for (i in which(is.na(h[,"x"]))) {
    if (i == 1L) {
      h$x[i] <- 0
      h$y[i] <- 0
    } else {
      h$x[i] <- h$x[i-1L]
      h$y[i] <- h$y[i-1L]
    }
  }
  h[h$t <= duration,,drop=FALSE]
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
  if (movement < 0)
    stop("movement can not be negative")
  if (any(mixture < 0))
    stop("mixture must not be negative")
  if (move_rate < 0)
    stop("move_rate must not be negative")
  if (vocal_rate < 0)
    stop("vocal_rate must not be negative")
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

get_events <- function(x, vocal_only=TRUE, tlim=NULL) {
  if (sum(x$abundance) == 0)
    return(data.frame(
      x=numeric(0),
      y=numeric(0),
      t=numeric(0),
      v=numeric(0),
      i=numeric(0)
    ))
  if (is.null(tlim))
    tlim <- c(0, x$duration)
  tlim <- pmin(pmax(0, tlim), x$duration)
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
  z <- z[z$t %[]% tlim,,drop=FALSE]
  z
}
points.bsims_events <-
function(x, vocal_only=TRUE, tlim=NULL, ...) {
  points(get_events(x, vocal_only, tlim), ...)
  invisible(x)
}
lines.bsims_events <-
function(x, tlim=NULL, ...) {
  if (is.null(tlim))
    tlim <- c(0, x$duration)
  tlim <- pmin(pmax(0, tlim), x$duration)
  N <- length(x$events)
  for (i in seq_len(N)) {
    xy <- cbind(
      x=x$nests$x[i] + x$events[[i]]$x,
      y=x$nests$y[i] + x$events[[i]]$y)
    xy <- xy[x$events[[i]]$t %[]% tlim,,drop=FALSE]
    lines(xy, ...)
  }
  invisible(x)
}
plot.bsims_events <-
function(x, tlim=NULL,
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
    lines(x, tlim=tlim, col=col_move, lty=lty_move, lwd=lwd_move, ...)
  if (!is.na(pch_vocal))
    points(x, vocal_only=TRUE, tlim=tlim,
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
  if (any(is.infinite(b)))
    stop("b must be finite")
  if (length(b) != length(tau)-1L)
    stop("length(b) must equal length(tau)-1")
  if (length(b) != length(unique(b)))
    stop("values in b must be unique")
  b <- sort(b)
  b[b==0] <- .Machine$double.eps
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
  vocal_only=TRUE, # should we detect visuals or just vocals?
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
    sobs <- cut(xy[1L], st, labels=FALSE)
  }
  for (i in seq_len(N)) {
    z <- x$events[[i]]
    ## bird position
    xxb <- x$nests$x[i] + z$x
    yyb <- x$nests$y[i] + z$y
    ## distance from observer
    xx <- xxb - xy[1L]
    yy <- yyb - xy[2L]
    z$d <- sqrt(xx^2 + yy^2)
    ## repel inds within repel distance (0 just for temp object z)
    z$v[z$d < repel] <- 0
    ## NA is placeholder for these vocalizations in return object
    x$events[[i]]$v[z$d < repel] <- NA
    if (vocal_only) {
      keep <- z$v > 0
      z <- z[keep,,drop=FALSE]
      xx <- xx[keep]
      yy <- yy[keep]
    }
    ## this is where HER attenuation somes in
    if (att) {
      theta <- atan2(yy, xx) # angle in rad
      sbrd <- cut(xxb, st, labels=FALSE)
      for (j in seq_len(nrow(z))) {
        ## order tau as HEREH
        TAU <- tau[c(1,2,3,2,1)][sobs:sbrd[j]]
        ## calculate distance breaks from x and theta
        if (length(TAU) == 1) {
          #b <- numeric(0)
          q <- dist_fun(z$d[j], TAU)
        } else {
          ## this gives breaks along x axis
          if (sobs < sbrd[j]) { # bird right of observer
            stj <- st[(sobs+1):sbrd[j]]
          } else { # bird left of observer
            stj <- st[sobs:(sbrd[j]+1)]
          }
          ## breaks as radial distance: r=x/cos(theta)
          b <- (stj - xy[1L]) / cos(theta[j])
          ## calculate q
          q <- dist_fun2(z$d[j], TAU, dist_fun, b)
        }
      }
    } else {
      q <- dist_fun(z$d, tau)
    }
    #u <- runif(length(z$d))
    #z$det <- ifelse(u <= q, 1, 0) # detected
    z$det <- rbinom(length(z$d), size=1, prob=q)
    z <- z[z$det > 0,,drop=FALSE]
    ## error is shown where detected, NA when not detected
    x$events[[i]]$d <- z$d[match(rownames(x$events[[i]]), rownames(z))]
  }
  x$xy <- xy
  x$tau <- tau
  x$repel <- repel
  x$vocal_only <- vocal_only
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
    "total duration: ", x$duration, "\n  detected: ", ndet,
    ifelse(x$vocal_only, " heard", " seen/heard"), "\n", sep="")
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
  ## angle in degrees counter clockwise from x axis right
  z$a <- 180 * atan2(z$x, z$y) / pi
  z$a[z$a < 0] <- 360+z$a[z$a < 0]
  ## observer position
  attr(z, "observer") <- x$xy
  z
}
points.bsims_detections <-
function(x, first_only=TRUE, ...) {
  points(get_detections(x, first_only)[,c("x", "y")], ...)
  invisible(x)
}
lines.bsims_detections <-
function(x, first_only=TRUE, tlim=NULL, ...) {
  if (is.null(tlim))
    tlim <- c(0, x$duration)
  z <- get_detections(x, first_only)
  z <- z[z$t %[]% tlim,,drop=FALSE]
  xy <- z[,c("x", "y")]
  segments(
    x0=rep(x$xy[1L], nrow(xy)),
    y0=rep(x$xy[2L], nrow(xy)),
    x1=xy[,1L], y1=xy[,2L], ...)
  invisible(x)
}

plot.bsims_detections <-
function(x, first_only=TRUE, tlim=NULL,
pch_nest=3, col_nest="darkgreen", cex_nest=1,
pch_vocal=21, col_vocal="blue", cex_vocal=0.5,
lty_move=1, col_move="orange", lwd_move=1,
lty_det=1, col_det="black", lwd_det=1,
...) {
  op <- par(xpd = TRUE)
  on.exit(par(op))
  xx <- x
  class(xx) <- c("bsim", "bsims_events")
  plot(xx, tlim=tlim,
    pch_nest=pch_nest, col_nest=col_nest, cex_nest=cex_nest,
    pch_vocal=pch_vocal, col_vocal=col_vocal, cex_vocal=cex_vocal,
    lty_move=lty_move, col_move=col_move, lwd_move=lwd_move, ...)
  if (!is.na(lty_det))
    lines(x, first_only, tlim=tlim, col=col_det, lty=lty_det, ...)
  invisible(x)
}


## lognormal parametrized as mean (ybar) and SDlog
rlnorm2 <- function(n, mean = exp(0.5), sdlog = 1) {
  rlnorm(n, log(mean) - sdlog^2/2, sdlog)
}
#summary(rlnorm2(10^6, 1.3, 0.5))

bsims_transcribe <- function(
  x,
  tint=NULL,
  rint=Inf,
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
  detall <- get_detections(x, first_only=FALSE)
  detall <- detall[detall$d <= max(rint),,drop=FALSE]
  if (error < 0)
    stop("error must be >= 0")
  derr <- if (error > 0)
    rlnorm2(nrow(detall), detall$d, error) else detall$d
  detall$error <- derr - detall$d
  rLAB <- paste0(c(0, round(100*rint[-length(rint)])),
    ifelse(is.finite(rint), paste0("-", round(100*rint)), "+"), "m")
  tLAB <- paste0(c(0, round(tint[-length(tint)], 2)), "-", tint, "min")
  detall$rint <- factor(rLAB[cut(derr, c(0, rint), labels=FALSE,
    include.lowest=TRUE)], rLAB)
  detall$tint <- factor(tLAB[cut(detall$t, c(0, tint), labels=FALSE,
    include.lowest=TRUE)], tLAB)

  ## count 1st detections over whole duration
  det <- detall
  if (first_only)
    det <- det[!duplicated(det$i),,drop=FALSE]
  xt <- as.matrix(Xtab(~ rint + tint, det))

  ## count 1st detections in visits (intervals)
  vis <- xt
  vis[] <- 0
  for (i in tLAB) {
    det2 <- detall[detall$tint == i,,drop=FALSE]
    if (first_only)
      det2 <- det2[!duplicated(det2$i),,drop=FALSE]
    vis <- vis + as.matrix(Xtab(~ rint + tint, det2))
  }

  x$detections <- detall
  x$removal <- xt
  x$visits <- vis
  x$tint <- tint
  x$rint <- rint
  x$first_only <- first_only
  x$error <- error
  class(x) <- c("bsim", "bsims_transcript", "bsims_detections")
  x
}
print.bsims_transcript <- function(x, ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  her <- paste0(
    ifelse(A[1] > 0, "H", ""),
    ifelse(A[2] > 0, "E", ""),
    ifelse(A[3] > 0, "R", ""), collapse="")
  ndet <- if (sum(x$abundance) == 0)
    0 else sum(sapply(x$events, function(z) any(!is.na(z$d))))
  cat("bSims transcript\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her,
    "\n  total abunance: ", sum(x$abundance),
    "\n  ", ifelse(length(x$mixture) > 0, "mixture with ", ""),
    "total duration: ", x$duration, "\n  detected: ", ndet,
    ifelse(x$vocal_only, " heard", " seen/heard"),
    "\n  ", ifelse(x$first_only, "1st", "all"),
    " inds. [", paste0(gsub("min", "", levels(x$det$tint)), collapse=", "),
    " min] [", paste0(gsub("m", "", levels(x$det$rint)), collapse=", "), " m]\n", sep="")
  invisible(x)
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
library(magrittr)

## check that abundance is right
l <- bsims_init(10)
summary(replicate(1000, sum(bsims_populate(l, 10)$abundance)) - 10^3)

## check vocal rates: no mixture
phi <- 0.5
br <- c(3, 5, 10)
#br <- 1:10
l <- bsims_init(10)
p <- bsims_populate(l, 1)
a <- bsims_animate(p, vocal_rate=phi)
o <- bsims_detect(a, tau=Inf) # detect all
d <- get_detections(o, first_only=TRUE)
i <- cut(d$t, c(0, br), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
D1 <- matrix(br, nrow=1)
(phihat <- exp(cmulti.fit(Y1, D1, type="rem")$coef))
plot(stepfun(d$t, (0:nrow(d))/nrow(d)), do.points=FALSE, xlim=c(0,10))
curve(1-exp(-phi*x), add=TRUE, col=2)
points(br, cumsum(table(i))/sum(table(i)), cex=2, col=4)
curve(1-exp(-phihat*x), add=TRUE, col=4)

## check vocal rates: finite mixture
phi <- c(10, 0.5)
mix <- c(0.2, 0.8)
#br <- c(3, 5, 10)
br <- 1:10
l <- bsims_init(10)
p <- bsims_populate(l, 10)
a <- bsims_animate(p, vocal_rate=phi, mixture=mix)
o <- bsims_detect(a, tau=Inf) # detect all
d <- get_detections(o, first_only=TRUE)
i <- cut(d$t, c(0, br), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
D1 <- matrix(br, nrow=1)
cf <- cmulti.fit(Y1, D1, type="mix")$coef # log.phi, logit.c
(phihat <- exp(cf[1]))
(mixhat <- c(1-plogis(cf[2]), plogis(cf[2])))


plot(stepfun(d$t, (0:nrow(d))/nrow(d)), do.points=FALSE, xlim=c(0,10))
curve(1-mix[2]*exp(-phi[2]*x), add=TRUE, col=2)
points(br, cumsum(table(i))/sum(table(i)), cex=2, col=4)
curve(1-mixhat[2]*exp(-phihat*x), add=TRUE, col=4)


## EDR

## check vocal rates: no mixture
phi <- 1
tau <- 0.8
br <- c(0.5, 1, 1.5, Inf)
l <- bsims_init(10)
Y1 <- D1 <- NULL
for (i in 1:20) {
  p <- bsims_populate(l, 10)
  a <- bsims_animate(p, vocal_rate=phi)
  o <- bsims_detect(a, tau=tau) # detect all
  d <- get_detections(o, first_only=TRUE)
  i <- cut(d$d, c(0, br), include.lowest = TRUE)
  Y1 <- rbind(Y1, matrix(as.numeric(table(i)), nrow=1))
  D1 <- rbind(D1, matrix(br, nrow=1))
}
Y1
D1
(tauhat <- exp(cmulti.fit(Y1, D1, type="dis")$coef))

plot(stepfun(d$t, (0:nrow(d))/nrow(d)), do.points=FALSE, xlim=c(0,10))
curve(1-exp(-phi*x), add=TRUE, col=2)
points(br, cumsum(table(i))/sum(table(i)), cex=2, col=4)
curve(1-exp(-phihat*x), add=TRUE, col=4)

## simple simulations
phi <- 1
tau <- 1
br <- c(0.5, 1, 1.5, Inf)
n <- 1000
x <- runif(n, -5, 5)
y <- runif(n, -5, 5)
d <- sqrt(x^2 + y^2)
p <- exp(-d^2/tau^2)
k <- rbinom(n, 1, p)
plot(x, y, asp=1, col="grey")
points(x[k>0], y[k>0], pch=19)
abline(h=0,v=0,lty=2)

i <- cut(d[k>0], c(0, br), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
D1 <- matrix(br, nrow=1)
(tauhat <- exp(cmulti.fit(Y1, D1, type="dis")$coef))

phi <- 0.25
dur <- 3
tau <- 1
br <- c(0.5, 1, 1.5, Inf)

testf <- function(phi, dur, tau=1, br=c(0.5, 1, 1.5, Inf), n=100) {
  res <- NULL
  for (iii in seq_len(n)) {
    l <- bsims_init(10)
    p <- bsims_populate(l, 10)
    a <- bsims_animate(p, vocal_rate=phi, duration=dur)
    o <- bsims_detect(a, tau=tau) # detect all

    x <- o$nests$x
    y <- o$nests$y
    d <- sqrt(x^2 + y^2)
    p <- exp(-d^2/tau^2)
    k <- rbinom(length(p), 1, p)
    i <- cut(d[k>0], c(0, br), include.lowest = TRUE)
    table(i)
    Y1 <- matrix(as.numeric(table(i)), nrow=1)
    D1 <- matrix(br, nrow=1)
    tauhat1 <- exp(cmulti.fit(Y1, D1, type="dis")$coef)

    z <- get_detections(o, first_only=TRUE)
    i <- cut(z$d, c(0, br), include.lowest = TRUE)
    table(i)
    Y1 <- matrix(as.numeric(table(i)), nrow=1)
    D1 <- matrix(br, nrow=1)
    tauhat2 <- exp(cmulti.fit(Y1, D1, type="dis")$coef)

    res <- rbind(res, c(nest=tauhat1, first=tauhat2))
  }
  list(phi=phi, dur=dur, tau=tau, br=br, res=res)
}

v <- expand.grid(phi=c(0.25, 0.5, 1), dur=c(3,5,10))

tt <- list()
for (j in 1:nrow(v)) {
  cat(j)
  tt[[j]] <- testf(phi=v$phi[j], dur=v$dur[j], n=10)
}

ttt <- t(sapply(tt, function(z) colMeans(z$res)))
cbind(v, ttt)
par(mfrow=c(3,3), mar=c(3,3,3,2))
for (j in 1:9) {
  boxplot(tt[[j]]$res/tau, main=paste(v$phi[j], v$dur[j]),
          ylim=c(0,1.3))
  abline(h=1,col=2)
}

## TODO:
## - use 1st vocal
## - use all vocals -> not very realistic in practice
## - use 1 randomly chosen vocal -> not very realistic in practice
## OK - estimate tau based on 3, 5, and 10 min duration

## both phi and tau
phi <- 0.5
tau <- 0.8
dur <- 10
rbr <- c(0.5, 1, 1.5, Inf)
tbr <- c(3, 5, 10)
l <- bsims_init(10)
p <- bsims_populate(l, 10)
a <- bsims_animate(p, vocal_rate=phi, duration=dur)
o <- bsims_detect(a, tau=tau) # detect all
x <- bsims_transcribe(o, tint=tbr, rint=rbr)
Y1 <- matrix(colSums(x$counts), nrow=1)
D1 <- matrix(x$tint, nrow=1)
Y2 <- matrix(rowSums(x$counts), nrow=1)
D2 <- matrix(x$rint, nrow=1)
exp(cmulti.fit(Y1, D1, type="rem")$coef)
exp(cmulti.fit(Y2, D2, type="dis")$coef)

## all vocals: OK
z <- get_detections(o, first_only=FALSE)
i <- cut(z$d, c(0, rbr), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
D1 <- matrix(rbr, nrow=1)
exp(cmulti.fit(Y1, D1, type="dis")$coef)
## random vocal: biased for large phi & dur
z <- z[sample(nrow(z)),]
z <- z[!duplicated(z$i),]
i <- cut(z$d, c(0, rbr), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
exp(cmulti.fit(Y1, D1, type="dis")$coef)
## 1st vocal: biased for large phi & dur
z <- get_detections(o, first_only=TRUE)
i <- cut(z$d, c(0, rbr), include.lowest = TRUE)
table(i)
Y1 <- matrix(as.numeric(table(i)), nrow=1)
exp(cmulti.fit(Y1, D1, type="dis")$coef)

## 3, 5, 10 min based EDR estimation: shorter the better
Y3 <- matrix(x$counts[,1], nrow=1)
Y5 <- matrix(rowSums(x$counts[,1:2]), nrow=1)
Y10 <- matrix(rowSums(x$counts), nrow=1)
D <- matrix(x$rint, nrow=1)
exp(cmulti.fit(Y3, D, type="dis")$coef)
exp(cmulti.fit(Y5, D, type="dis")$coef)
exp(cmulti.fit(Y10, D, type="dis")$coef)


## avoid R stratum
l <- bsims_init(10, 0.5, 0.5)
plot(l)
p <- bsims_populate(l, c(1,1,0))
plot(p)
a <- bsims_animate(p, movement=0.25, move_rate=1, avoid="R")
plot(a, tlim=c(5,10))
o <- bsims_detect(a, c(0.5,0))
plot(o, tlim=c(5,10))
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
plot(o, pch_vocal=NA)

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


col2hex <- function(col, alpha = FALSE) {
  rgb <- col2rgb(col, alpha)
  if (alpha) {
    apply(rgb, 2, function(z) rgb(z[1], z[2], z[3], z[4], maxColorValue=255))
  } else {
    apply(rgb, 2, function(z) rgb(z[1], z[2], z[3], maxColorValue=255))
  }
}
#col2hex(c(blu = "royalblue", reddish = "tomato"), FALSE)
#col2hex(c(blu = "royalblue", reddish = "tomato"), TRUE)
color_fade <- function(col, n) {
  rgb <- col2rgb(col[1L], alpha=FALSE)
  sapply(seq(0, 255, length.out = n), function(z)
    rgb(rgb[1], rgb[2], rgb[3], z, maxColorValue=255))
}
#plot(1:10, col=color_fade("red", 10), pch=19)


library(magick)

set.seed(1234)
l <- bsims_init(15, 0.1, 0.5, 2.5)
#l <- bsims_init(10, 0.1, 0.5)
p <- bsims_populate(l, c(2, 1.5, 0))
a <- bsims_animate(p, duration=60, movement=0.2, move_rate=2, vocal_rate=2, avoid="R")
o <- bsims_detect(a, xy=c(2.5, 1), tau=c(1:3)) # detect all
plot(o, pch_nest=NA, pch_vocal=NA, first_only=FALSE, tlim=c(0,60))
mtext("bSims: highly scientific and utterly addictive", 1, -2)

plot(o, pch_nest=NA, pch_vocal=NA, first_only=FALSE, tlim=c(0,60), xlim=c(-4,-2), ylim=c(-4, -2))
plot(o, pch_nest=NA, pch_vocal=NA, first_only=FALSE, tlim=c(0,60), xlim=c(-4,4), ylim=c(-4, 4))

chr0 <- function(x, n) {
  x <- as.character(x)
  paste0(paste0(rep(0, n - nchar(x)), collapse=""), x, collapse="")
}
## quadratic ease in-out: https://gist.github.com/gre/1650294
g <- function(from, to, length.out) {
  t <- seq(0, 1, length.out = length.out)
  a <- ifelse(t < 0.5, 2*t*t, -1+(4-2*t)*t)
  a * (to-from) + from
}
fps <- 19
lgt <- 6
ti <- 3
nstep <- lgt*fps

x0 <- g(-4, -4, length.out = nstep)
x1 <- g(-2, 4, length.out = nstep)
y0 <- g(-4, -4, length.out = nstep)
y1 <- g(-2, 4, length.out = nstep)
t0 <- g(0, o$duration-ti, length.out = nstep)
t1 <- g(ti, o$duration, length.out = nstep)

par(mar=rep(0,4))
for (i in 1:nstep) {
  png(paste0("temp/plot", chr0(i,3), ".png"), width=400, heigh=400)
  plot(o, pch_nest=NA, pch_vocal=NA, first_only=FALSE,
    tlim=c(t0[i],t1[i]), xlim=c(x0[i],x1[i]), ylim=c(y0[i],y1[i]))
  mtext("bSims: highly scientific and utterly addictive", 1, 0)
  dev.off()
}

im <- image_read(paste0("temp/", list.files("temp", pattern=".png")))
an <- image_animate(im)
image_write(an, "temp/bsims.gif")
