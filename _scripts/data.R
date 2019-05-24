library(mefa4)
load("~/Dropbox/courses/aos-2019-anchorage/data/ab-birds-all-2018-11-29.RData")

## find center point not that is in grassland
zz <- droplevels(dd[dd$CMETHOD=="RF" & substr(as.character(dd$SS),1,2) != "OG",])
tmp <- strsplit(as.character(zz$SS), "_")
zz$ABMIsite <- sapply(tmp, "[[", 1)
zz$ABMIbirdpt <- sapply(tmp, "[[", 2)
n <- table(zz$ABMIsite)
n <- n[n==9]
zz$All9 <- zz$ABMIsite %in% names(n)
nn <- sum_by(zz$NRNAME == "Grassland" | zz$Y < 50, zz$ABMIsite)
zz$NotGr <- zz$ABMIsite %in% rownames(nn)[nn[,"x"] == 0]

table(ngr=zz$NotGr, a=zz$All9)/9

nam <- unique(zz[zz$All9, "ABMIsite"])

dd$ABMIsite <- zz$ABMIsite[match(rownames(dd), rownames(zz))]
dd$ABMIsite[is.na(dd$ABMIsite)] <- ""

rn <- rownames(dd[dd$ABMIsite %in% nam,])

dat <- dd[rn,]
veg <- vc1[rn,]
veg <- veg / rowSums(veg)
dat$pr <- as.numeric(veg[,"RoadHardSurface"])
dat$ph <- as.numeric(veg[,"RoadVegetatedVerge"])
## there are roads without verge
with(dat[dat$pr > 0,], plot(pr, pr+ph))
abline(0,1)
## need to have broad forest types identified (decid+mixed, conif+pine)
## and match road PC with off-road PCs at same site
## pr * ph * veg + (1|site)
## check

dat$pdec <- rowSums(veg[,grepl("Decid", colnames(veg)) | grepl("Mixedwood", colnames(veg))])
dat$pcon <- rowSums(veg[,grepl("Spruce", colnames(veg)) |
    grepl("Pine", colnames(veg))])
dat$plocon <- rowSums(veg[,grepl("TreedFen", colnames(veg)) |
    grepl("TreedBog", colnames(veg))])

with(dat[dat$pr > 0,], plot(pdec, pcon))
with(dat, plot(pdec, pcon))
with(dat[dat$pr > 0,], points(pdec, pcon, col=2))
dat$roadside <- dat$pr > 0.01
dat$offroad <- dat$pr == 0
dat$decid <- dat$pdec > 0.75
dat$upfor <- dat$pdec + dat$pcon
dat$pfor <- dat$pdec + dat$pcon + dat$plocon
dat$popen <- rowSums(veg[,c("GrassHerb", "Shrub")])
dat$open <- dat$popen > 0.75
dat <- dat[!is.na(dat$roadside),]
dat$site_year <- paste0(substr(as.character(dat$SS), 1, nchar(as.character(dat$SS))-2), "_",
                        dat$YEAR)


with(dat, plot(open, pfor))

dat <- dat[dat$roadside | dat$offroad,]
dat <- dat[dat$decid | dat$open,]
#dat <- dat[dat$for,]
table(table(dat$site_year))

dat$ok <- FALSE
for (i in unique(dat$site_year)) {
  tmp <- dat[dat$site_year == i,]
  if (nrow(tmp) > 1) {
    if (any(tmp$roadside) && any(tmp$offroad))
      dat$ok[dat$site_year == i] <- TRUE
  }
}
dat <- dat[dat$ok,]
dat$rrr <- ifelse(dat$roadside, 1, 0)
dat$fff <- ifelse(dat$pfor, 1, 0)
table(dat$site_year,dat$roadside)

Y <- yy[rownames(dat),]
Y <- Y[,colSums(Y) > 0]


library(lme4)
library(mgcv)

spp <- "TEWA"
dat$count <- as.numeric(Y[,spp])

m <- gam(count ~ s(pr) + s(ph), data=dat, family=poisson)
plot(m)
m <- step(glm(count ~ rrr * fff, data=dat, family=poisson))
m <- glm(count ~ rrr + fff, data=dat, family=poisson)
#m <- glmer(count ~ rrr + (1|site_year), data=dat, family=poisson)
summary(m)

