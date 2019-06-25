library(intrval)
x <- read.csv("_data/bayne-2016/condor-suppl.csv", stringsAsFactors=FALSE)
x1 <- x[,1:4]
x1[] <- lapply(x1, as.factor)
rownames(x1) <- x1$AOUcode
x2 <- x[,5:13]
x3 <- x[,14:22]

f1 <- function(i) {
    z <- t(sapply(strsplit(x2[,i], "/"), as.integer))
    colnames(z) <- paste0(colnames(x2)[i], c("_Impact", "_Control"))
    z
}
x2 <- do.call(cbind, lapply(1:ncol(x2), f1))
rownames(x2) <- x1$AOUcode

f21 <- function(z) {
    if (z == "")
        return(c(NA, NA, NA))
    zz <- strsplit(z, " ")
    zzz <- strsplit(zz[[1]][2], "-")
    as.numeric(gsub("[^0-9\\.]", "", c(zz[[1]][1], zzz[[1]])))
}
f2 <- function(i) {
    z <- t(sapply(x3[,i], f21))
    colnames(z) <- paste0(colnames(x3)[i], c("", "_Lo", "_Hi"))
    z
}
x3 <- do.call(cbind, lapply(1:ncol(x3), f2))
rownames(x3) <- x1$AOUcode
xx <- data.frame(x1, x2, x3)

spp <- "ALFL"

m <- matrix(x3[spp,], nrow=3)
m0 <- matrix(m[1,], 3, 3)
mL <- matrix(m[2,], 3, 3)
mU <- matrix(m[3,], 3, 3)
mSD <- (log(mU) - log(mL)) / qnorm(1-0.1/2)
dimnames(m0) <- dimnames(mL) <- dimnames(mU) <-
    dimnames(mSD) <- list(c(50, 100, Inf), c("S", "P", "W"))

i <- "50"
j <- "S"

vals <- rnorm(10^6, log(m0[i,j]), mSD[i,j])
plot(density(exp(vals)))

col <- c('#66c2a5','#fc8d62','#e78ac3')

pdf("spp.pdf", onefile=TRUE, height=5, width=7)
for (spp in rownames(x1)) {
    m <- matrix(x3[spp,], nrow=3)
    m0 <- matrix(m[1,], 3, 3)
    mL <- matrix(m[2,], 3, 3)
    mU <- matrix(m[3,], 3, 3)
    main <- paste(x1[spp, "CommonName"], x1[spp, "Habitat"])
    colnames(m0) <- c("Seismic", "Pipeline", "Wellpad")
    m0[is.na(m0)] <- 0
    mL[is.na(mL)] <- 0
    mU[is.na(mU)] <- 0
    v <- barplot(m0,beside=TRUE, main=main, col=col, ylim=c(0, min(10, max(mU))))
    abline(h=1, col='#8da0cb')
    for (k in 1:9) {
        int <- c(mL[k], mU[k])
        lwd <- if (1 %[]% int) 1 else 3
        lines(c(v[k], v[k]), int, lwd=lwd)
    }
}
dev.off()

x3x <- x3[,!grepl("_", colnames(x3))]

cc <- c(1,2,3,NA,4,5,6,NA,7,8,9)

pdf("sppbox.pdf", onefile=TRUE, height=5, width=7)
boxplot(x3x[,cc], ylim=c(0, 3), main="", col=c(col,NA), axes=FALSE)
axis(2)
axis(1,c(2,6,10),c("Seismic", "Pipeline", "Wellpad"), tick=F)
abline(h=1, col='#8da0cb')
boxplot(x3x[x1$Habitat!="OLD",cc], ylim=c(0, 3), main="Not OLD", col=c(col,NA), axes=FALSE)
axis(2)
axis(1,c(2,6,10),c("Seismic", "Pipeline", "Wellpad"), tick=F)
abline(h=1, col='#8da0cb')
for (i in c("OLD", "ALL", "SHR")) {
    boxplot(x3x[x1$Habitat==i,cc], ylim=c(0, 3), main=i, col=c(col,NA), axes=FALSE,
        ylab="Relative effect")
    axis(2)
    axis(1,c(2,6,10),c("Seismic", "Pipeline", "Wellpad"), tick=F)
    abline(h=1, col='#8da0cb')
    legend("topleft", fill=col, legend=c("50m", "100m", "Unlimited"), bty="n", title="Truncation distance")
}
dev.off()

op <- par(mfrow=c(2,1))
i <- "OLD"
boxplot(x3x[x1$Habitat==i,cc], ylim=c(0, 3), main="Old forest species",
    col=c(col,NA), axes=FALSE, ylab="Relative effect")
axis(2)
axis(1,c(2,6,10),c("Seismic", "Pipeline", "Wellpad"), tick=F)
abline(h=1, col='#8da0cb')
legend("topleft", fill=col, legend=c("50m", "100m", "Unlimited"), bty="n", title="Truncation distance")

i <- "SHR"
boxplot(x3x[x1$Habitat==i,cc], ylim=c(0, 3), main="Shrub species",
    col=c(col,NA), axes=FALSE, ylab="Relative effect")
axis(2)
axis(1,c(2,6,10),c("Seismic", "Pipeline", "Wellpad"), tick=F)
abline(h=1, col='#8da0cb')
legend("topleft", fill=col, legend=c("50m", "100m", "Unlimited"), bty="n", title="Truncation distance")
par(op)


