#' # bSims: bird simulations
#'
#' ## Landscape initialization
#'
library(intrval)
library(MASS)
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
  col=c("green", "yellow", "grey"), ...) {
  A <- diff(x$strata) * diff(range(x$strata))
  A <- c(h=A[1]+A[5], e=A[2]+A[4], r=A[3])
  names(A) <- c("H", "E", "R")
  plot(0, type="n", xlim=range(x$box[,"x"]), ylim=range(x$box[,"y"]),
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
  ...) {
  if (!inherits(x, "bsims_landscape"))
    stop("x must be a bsims_landscape object")
  A <- diff(x$strata) * diff(range(x$strata))
  D <- rep(density, 3)[c(1,2,3,2,1)]
  lambda <- A * D
  if (is.null(abund_fun))
    abund_fun <- function(lambda, ...) rpois(1, lambda)
  N <- sapply(lambda, abund_fun, ...)
  d <- NULL
  for (i in 1:5) {
    if (xy_process == "poisson") {
      z <-data.frame(
        acceptreject(n=N[i], f=xy_fun,
          x0=x$strata[i], x1=x$strata[i+1],
          y0=x$strata[1], y1=x$strata[6],
          m=margin,
          maxit=maxit,
          fail=fail),
        s=rep(i, N[i]))
    }
    # add here spatial non-randomness
    d <- rbind(d, z)
  }
  d$i <- 1:nrow(d)
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
function(x, add=FALSE, ...) {
  if (!add)
    plot(0, type="n", xlim=range(x$box[,"x"]), ylim=range(x$box[,"y"]),
      xlab="", ylab="", axes=FALSE, asp=1, ...)
  points(x, ...)
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
events <-
function(vocal_rate=1, move_rate=1, duration=10, movement=0) {
  ev <- cumsum(timetoevent(vocal_rate, duration))
  em <- cumsum(timetoevent(move_rate, duration))
  iv <- rep(1, length(ev))
  im <- rep(0, length(em))
  dv <- matrix(NA, length(ev), 2)
  dm <- MASS::mvrnorm(length(em), c(0, 0), diag(movement^2, 2, 2))
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

bsims_animate <- function(
  x, # population object
  vocal_rate=1, # phi /min
  move_rate=1, #movement
  duration=10,
  movement=0, # SD for 2D kernel
  mixture=1, # finite mixture group proportions
  ...) {
  if (!inherits(x, "bsims_population"))
    stop("x must be a bsims_population object")
  if (any(mixture < 0))
    stop("mixture must not be negative")
  K <- length(mixture)
  G <- paste0("G", 1:K)
  P <- mixture / sum(mixture)
  vocal_rate <- rep(vocal_rate, K)[1:K]
  move_rate <- rep(move_rate, K)[1:K]
  names(P) <- names(vocal_rate) <- names(move_rate) <- G
  N <- sum(x$abundance)
  g <- sample(G, N, replace=TRUE, prob=P)
  x$nests$g <- factor(g, G)
  Events <- list()
  for (i in seq_len(N)) {
    Events[[i]] <- events(
      vocal_rate=vocal_rate[g[i]],
      move_rate=move_rate[g[i]],
      duration=duration, movement=movement)
  }
  x$vocal_rate <- vocal_rate
  x$move_rate <- move_rate
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
  cat("bSims population\n  ",
    round(x$extent/10, 1), " km x ", round(x$extent/10, 1),
    " km\n  stratification: ", her,
    "\n  total abunance: ", sum(x$abundance),
    "\n  ", ifelse(length(x$mixture) > 0, "mixture with ", ""),
    "total duration: ", x$duration, "\n", sep="")
  invisible(x)
}

points.bsims_events <-
function(x, vocal_only=TRUE, ...) {
  N <- length(x$events)
  for (i in seq_len(N)) {
    xy <- cbind(
      x=x$nests$x[i] + x$events[[i]]$x,
      y=x$nests$y[i] + x$events[[i]]$y)
    if (vocal_only)
      xy <- xy[x$events[[i]]$v > 0,,drop=FALSE]
    points(xy, ...)
  }
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
function(x, add=FALSE,
col.line=1, col.point=1,
lty=1, pch=21, ...) {
  if (!add)
    plot(0, type="n", xlim=range(x$box[,"x"]), ylim=range(x$box[,"y"]),
      xlab="", ylab="", axes=FALSE, asp=1, ...)
  N <- length(x$events)
  for (i in seq_len(N)) {
    xy <- cbind(
      x=x$nests$x[i] + x$events[[i]]$x,
      y=x$nests$y[i] + x$events[[i]]$y)
    if (!is.na(lty))
      lines(xy, ...)
    if (!is.na(pch))
      points(xy[x$events[[i]]$v > 0,,drop=FALSE], pch=pch, ...)
  }
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


l <- bsims_init(3)
p <- bsims_populate(l, c(2,1,0))
x <- bsims_animate(p, movement=0.1)
plot(l)
lines(x, col="grey")
points(x, col=2, cex=0.4)
points(p, pch=3, cex=0.6)


rr <- 1
tt <- timetoevent(rr, 10)
op <- par(mfrow=c(1,2))
plot(ecdf(tt))
curve(1-exp(-rr*x), add=TRUE, col=2)

plot(stepfun(sort(tt), 0:length(tt)/length(tt)))
curve(1-exp(-rr*x), add=TRUE, col=2)
par(op)

