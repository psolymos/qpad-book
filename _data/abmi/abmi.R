#' ---
#' title: "ABMI bird data set processing"
#' output: pdf_document
#' ---
#'
library(mefa4)
#knitr::opts_chunk$set(eval=FALSE)
if (interactive()) {
  od <- setwd("_data/abmi")
}
#' Read in table
d <- read.csv("abmi-song-meter-2015-2017.csv")
str(d)
#' Format start date and time
tmp <- paste(d$Recording_Date, d$Recording_Time)
d$Start <- strptime(tmp, "%d-%b-%y %H:%M:%S")
#' Add total duration based on listening methodology
d$Duration <- NA
d$Duration[d$Method %in% c("11", "14")] <- 3
d$Duration[d$Method %in% c("12", "13")] <- 1
#' Define primary key as site+year+quadrant
d$pkey <- as.factor(paste0(d$ABMI_Site, "_", d$Year, "_", d$Quadrant))
#' Visit add date and time to pkey
d$visit <- as.factor(paste0(d$pkey, "_", gsub(" ", "_", as.character(d$Start))))
#' Use CamelCase species IDs
d$SpeciesID <- d$Common_Name
levels(d$SpeciesID) <- nameAlnum(levels(d$SpeciesID), capitalize="mixed", collapse="")
d$SpeciesID <- droplevels(d$SpeciesID)
#' First detection interval
d$int1 <- ifelse(d$Interval_1_1_minute == "VNA", NA, as.integer(d$Interval_1_1_minute))
d$int2 <- ifelse(d$Interval_2_1_minute == "VNA", NA, as.integer(d$Interval_2_1_minute))
d$int3 <- ifelse(d$Interval_3_1_minute == "VNA", NA, as.integer(d$Interval_3_1_minute))
tmp <- col(d[,c("int1", "int2", "int3")])
tmp[is.na(d[,c("int1", "int2", "int3")])] <- Inf
tmp2 <- find_min(tmp)
tmp2$value[is.infinite(tmp2$value)] <- NA
d$det1 <- factor(c("0-1min", "1-2min", "2-3min")[tmp2$value],
  c("0-1min", "1-2min", "2-3min"))
#' Time of year and time of day variables
d$ToY <- d$Start$yday
d$ToYc <- as.integer(cut(d$ToY, c(0, 105, 120, 140, 150, 160, 170, 180, 365)))
d$ToD <- d$Start$hour + d$Start$min / 60
d$ToDx <- round(d$ToD, 0)
d$ToDc <- as.factor(ifelse(d$ToDx < 4, "Midnight", "Morning"))
#' Save object
keep <- !is.na(d$Duration) & !is.na(d$Start)
abmi <- droplevels(d[keep,])
if (interactive()) {
  save(abmi, file="abmi.rda")
  setwd(od)
}
#' Make survey x interval tables for species based on 3-min samples:
#' need to drop visits that are <3 min
y3 <- Xtab(~ visit + det1 + SpeciesID,
  data=abmi[abmi$Duration == 3 & !(abmi$SpeciesID %in% c("NONE","SNI", "VNA", "DNC", "PNA")),],
  drop.unused.levels=TRUE)
x3 <- droplevels(nonDuplicated(abmi, visit, TRUE)[rownames(y3[[1]]),
  c("pkey", "visit", "ToY", "ToYc", "ToD", "ToDx", "ToDc")])

library(detect)
ToY <- seq(min(x1$ToY), max(x1$ToY), 1)

spp <- "Ovenbird"
#spp <- "BorealChickadee"

D <- as.matrix(y3[[1]])
D[,1] <- 1
D[,2] <- 2
D[,3] <- 3
Y <- as.matrix(y3[[spp]])
Y <- Y[x3$ToDc == "Morning",]
D <- D[x3$ToDc == "Morning",]
z <- cmulti(Y|D ~ ToY+I(ToY^2)+I(ToY^3),
  x3[x3$ToDc == "Morning",], type="rem")
X <- model.matrix(~ ToY+I(ToY^2)+I(ToY^3),
  data.frame(ToY=seq(min(x1$ToY), max(x1$ToY), 1)))
phi <- exp(X %*% coef(z))
p0 <- 1-exp(-1*phi)

z <- cmulti(Y|D ~ ToY+I(ToY^2)+I(ToY^3),
  x3[x3$ToDc == "Morning",], type="mix")
summary(z)
chat <- plogis(X %*% coef(z)[-1])
phi <- exp(coef(z)[1])
p1 <- 1-chat*exp(-1*phi)
plot(p1 ~ ToY, type="l", ylim=c(0,1))
lines(p0 ~ ToY, col=4)

#' Make survey x species table based on 1st 1-min:
#' exclude counts from >1min but don't drop visit labels
y1 <- Xtab(~ visit + SpeciesID,
  data=abmi[is.na(abmi$det1) | abmi$det1 == "0-1min",],
  cdrop=c("NONE","SNI", "VNA", "DNC", "PNA"))
x1 <- droplevels(nonDuplicated(abmi, visit, TRUE)[rownames(y1),
  c("pkey", "visit", "ToY", "ToYc", "ToD", "ToDx", "ToDc")])

x1$y <- ifelse(y1[,spp] > 0, 1, 0)
tmp <- with(x1[x1$ToY >= 150, ], table(pkey, y))
tmp <- tmp[tmp[,"1"] > 0,]
x1$occ <- x1$pkey %in% rownames(tmp)

m <- mgcv::gam(y ~ s(ToY), x1[x1$ToDc == "Morning" & x1$occ,], family=binomial)
p <- predict(m, newdata=data.frame(ToY=ToY), type="response")

m <- mgcv::gam(y ~ s(ToY), x1[x1$ToDc == "Morning" & x1$occ &
    x1$visit %in% levels(x3$visit),], family=binomial)
p3 <- predict(m, newdata=data.frame(ToY=ToY), type="response")


#' Multi species
s <- colSums(y1[x1$ToY >= 150, ] > 0)
SPP <- colnames(y1)[s >= 100]
ToY <- seq(min(x1$ToY), max(x1$ToY), 1)
P <- matrix(0, length(ToY), length(SPP))
colnames(P) <- SPP

for (spp in SPP) {
  x1$y <- ifelse(y1[,spp] > 0, 1, 0)
  tmp <- with(x1[x1$ToY >= 150, ], table(pkey, y))
  tmp <- tmp[tmp[,"1"] > 0,]
  x1$occ <- x1$pkey %in% rownames(tmp)
  m <- mgcv::gam(y ~ s(ToY), x1[x1$ToDc == "Morning" & x1$occ,], family=poisson)
  P[,spp] <- predict(m, newdata=data.frame(ToY=ToY), type="response")
}

matplot(ToY, P, type="l", ylim=c(0,1), lty=1, col="#00000044")
rug(jitter(x1$ToY), side=1, col=2)
rug(jitter(x3$ToY), side=3, col=4)
quantile(x1$ToY, c(0.025, 0.975))
quantile(x3$ToY, c(0.025, 0.975))

#' Removal model
library(survival)

x1$y <- ifelse(y1[,spp] > 0, 1, 0)
tmp <- with(x1[x1$ToY >= 150, ], table(pkey, y))
tmp <- tmp[tmp[,"1"] > 0,]
x1$occ <- x1$pkey %in% rownames(tmp)
x1 <- droplevels(x1[x1$occ,])
xx <- droplevels(abmi[abmi$SpeciesID == spp & !is.na(abmi$int1),])
#' Events are coded into a survival object:
#' we treat nondetections as censored events at 60 sec,
x1$time <- xx$int1[match(x1$visit, xx$visit)]
x1$time[is.na(x1$time)] <- 60
#' time cannot be 0, so we use 1 sec instead
x1$time[x1$time == 0] <- 1
#' We give time in minutes, so we get rate as events/min.
x1$sv <- Surv(x1$time/60, x1$y)
#' Fit a series of survival models
mods <- list(m0 = survreg(sv ~ 1, x1, dist="exponential"),
    m1 = survreg(sv ~ ToY, x1, dist="exponential"),
    m2 = survreg(sv ~ ToY+I(ToY^2), x1, dist="exponential"),
    m3 = survreg(sv ~ ToY+I(ToY^2)+I(ToY^3), x1, dist="exponential"))

aic <- data.frame(df=sapply(mods, function(z) length(coef(z))), AIC=sapply(mods, AIC))
aic$dAIC <- aic$AIC - min(aic$AIC)
aic
#' `survreg` fits accelerated failure models, not proportional
#' hazards models, so the coefficients are logarithms of ratios of
#' survival times, and a positive coefficient means longer survival.
mb <- mods[[which.min(aic$AIC)]]
summary(mb)
#' Survival times
summary(predict(mb, data.frame(ToY=ToY)))
#' Event rate per unit (1 min) time
summary(1/predict(mb, data.frame(ToY=ToY)))
#' Probability of at least 1 event per 10 time units (mins)
p5 <- 1-exp(-1*(1/predict(mb, data.frame(ToY=ToY))))
summary(p5)


plot(p ~ ToY, type="l", ylim=c(0,1))
lines(p3 ~ ToY, col=3)
lines(p1 ~ ToY, col=2)
lines(p0 ~ ToY, col=4)
lines(p5 ~ ToY, col=5)
rug(jitter(x1$ToY), side=1, col=2)
rug(jitter(x3$ToY), side=3, col=4)

## basic survival model simulation

phi <- 2 # rate
x <- rexp(1000, 2)
summary(x)
TT <- 1
io <- x > TT
x[x > TT] <- TT
# time: right censored data, this is the follow up time (time of 1st det)
# event: status indicator, normally 0=alive, 1=dead (dead=heard)
sv <- Surv(x, ifelse(io, 0, 1))
m <- survreg(sv ~ 1, dist="exponential")

summary(m)
exp(-coef(m))
summary(1/predict(m))
