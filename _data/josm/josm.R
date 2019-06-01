#' ---
#' title: "JOSM bird data set processing"
#' output: pdf_document
#' ---
#'
#' # Preamble
#'
library(mefa4)
#knitr::opts_chunk$set(eval=FALSE)
if (interactive())
  setwd("_data/josm")
#'
#' # Bird species
#'
#' From data set 1
s1 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2011-13_BirdSpecies-v1.csv",
  skip=17, encoding="latin1")
str(s1)
s1$ScientificName <- s1$Scientific.Name
#' From data set 2
s2 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2014_BirdSpecies-v1.csv",
  skip=18)
str(s2)
#' Bind the 2 tables and select columns of interest
s12 <- rbind(
  s1[,c("SpeciesID", "SpeciesName", "ScientificName")],
  s2[,c("SpeciesID", "SpeciesName", "ScientificName")]
)
#' Keep unique AOU codes and use `SpeciesID`s as row names
s12 <- s12[!duplicated(s12$SpeciesID),]
rownames(s12) <- s12$SpeciesID
summary(s12)
#'
#' # Sites
#'
#' From data set 1
x1 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2011-13_Sites-v1.csv",
  skip=17)
str(x1)
#' From data set 2
x2 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2014_Sites-v1.csv",
  skip=18)
str(x2)
#' Bind the 2 tables and select columns of interest, use `SiteID`s as row names
cn <- c("SiteID", "SurveyArea", "Longitude", "Latitude")
x12 <- rbind(x1[,cn], x2[,cn])
rownames(x12) <- x12$SiteID
summary(x12)
plot(Latitude ~ Longitude, x12, pch=".")
#'
#' # Surveys
#'
#' Stations are 300m apart at the site location.
#'
#' From data set 1
z1 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2011-13_SiteSurveyInfo-v1.csv",
  skip=17)
str(z1)
z1$Date <- as.Date(as.character(z1$Date), "%d/%m/%Y")
#' From data set 2
z2 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2014_SiteSurveyInfo-v1.csv",
  skip=18)
str(z2)
z2$Date <- as.Date(as.character(z2$Date), "%m/%d/%Y")
colnames(z2)[colnames(z2) == "PCSiteID"] <- "SiteID"
colnames(z2)[colnames(z2) == "O.N.Rain"] <- "OvernightRain"
#' Bind the 2 tables and select columns of interest, use `StationID`s as row names
cn <- c("Date", "SiteID", "StationID", "ObserverID", "TimeStart",
  "VisitID", "WindStart", "PrecipStart", "TempStart", "CloudStart",
  "WindEnd", "PrecipEnd", "TempEnd", "CloudEnd", "TimeFin", "Noise",
  "OvernightRain")
z12 <- rbind(z1[,cn], z2[,cn])
rownames(z12) <- z12$SiteID
levels(z12$OvernightRain) <-startsWith(toupper(levels(z12$OvernightRain)), "Y")
summary(z12)
#' Check if all `SiteID`s are accounted for
compare_sets(x12$SiteID,z12$SiteID)
#' Deal with date/time (some days and months are )
tmp <- paste(z12$Date, z12$TimeStart)
d1 <- as.POSIXlt(tmp, tz="America/Edmonton", format="%Y-%m-%d %H:%M:%OS AM")
tmp <- paste(z12$Date, z12$TimeFin)
d2 <- as.POSIXlt(tmp, tz="America/Edmonton", format="%Y-%m-%d %H:%M:%OS AM")
#' There are `NA` in `d1` so calculate start time from `d2` instead
summary(as.numeric(d1-d2))
z12$DateTime <- d2 - 600
#'
#' Calculate surise time
library(maptools)
srise <- sunriset(
  as.matrix(x12[match(z12$SiteID, x12$SiteID), c("Longitude", "Latitude")]),
  z12$DateTime,
  direction="sunrise",
  POSIXct.out=TRUE)
z12$SunRiseTime <- srise$time
z12$SunRiseFrac <- srise$day_frac
z12$TSSR <- as.numeric(z12$DateTime - z12$SunRiseTime) / (60*60*24)
#' Ordinal day
z12$OrdinalDay <- as.POSIXlt(z12$DateTime)$yday
z12$DAY <- z12$OrdinalDay / 365
plot(TSSR ~ jitter(DAY), z12, pch=".")
#'
#' # Species counts
#'
y1 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2011-13_SpeciesObservationLog-v1.csv",
  skip=17)
str(y1)
#' From data set 2
y2 <- read.csv(
  "Cause-Effect-Biodiversity-Monitoring-Landbirds-2014_SpeciesObservationLog-v1.csv",
  skip=18)
str(y2)
colnames(y2)[colnames(y2) == "X"] <- "ObservationID"
colnames(y2)[colnames(y2) == "PCSiteID"] <- "SiteID"
cn <- c("ObservationID", "SiteID", "StationID", "Species",
  "TimeInterval", "Direction", "Distance",
  "DetectType1", "DetectType2", "DetectType3",
  "Sex", "Age", "Activity1", "Activity2", "Activity3", "ActivityNote")
y12 <- rbind(y1[,cn], y2[,cn])
rownames(y12) <- y12$ObservationID
summary(y12)
for (i in colnames(y12))
  if (is.factor(y12[[i]])) {
    levels(y12[[i]]) <- toupper(levels(y12[[i]]))
    y12[[i]][y12[[i]] == ""] <- NA
    y12[[i]] <- droplevels(y12[[i]])
  }
#' Some tweaks to species labels
levels(y12$Species)[levels(y12$Species) == "YEWA"] <- "YWAR"
y12 <- y12[y12$Species %in% s12$SpeciesID,]
y12$Species <- droplevels(y12$Species)
#' Drop rows where some key variables are `NA`
y12 <- y12[!is.na(y12$TimeInterval) &
  !is.na(y12$Distance) &
  y12$Distance != 0 & # unknown distance
  !is.na(y12$DetectType1), ]
table(y12$Distance,y12$TimeInterval,useNA = "a")
#' Specify human readable time and distance intervals
y12$Dur <- as.factor(y12$TimeInterval)
levels(y12$Dur) <- c("0-3min", "3-5min", "5-10min")
y12$Dis <- as.factor(y12$Distance)
levels(y12$Dis) <-c("0-50m", "50-100min", "100+m", "100+m")
table(y12$Dur, y12$Dis, useNA = "a")
#' Some other columns of interest
table(y12$Sex, y12$Age, useNA = "a")
table(y12$DetectType1, useNA="a")
#' Check if all `SiteID`s and `StationID`s are accounted for
compare_sets(x12$SiteID, y12$SiteID)
compare_sets(z12$StationID, y12$StationID)
#'
#' # Spatial covariates
#'
#' These 1km resolution rasters follow Solymos et al. in review (PIF/PIX manuscript)
library(raster)
library(sp)
rr <- stack("landcover-hfi2016.grd")
#' Define CRS NAD83 for our sites
xy <- x12[,c("Longitude", "Latitude")]
coordinates(xy) <- ~ Longitude + Latitude
proj4string(xy) <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
xy <- spTransform(xy, proj4string(rr))
#' Plot the sites over the water layer
col <- colorRampPalette(c("lightgrey", "blue"))(100)
plot(rr[["Water"]], col=col, axes=FALSE, box=FALSE)
plot(xy, add=TRUE, pch=19, cex=0.5)
#'
#' Extract values at each point
v <- extract(rr, xy)
#'
#' # Assamble and save the output object
#'
josm <- list(
  sites=data.frame(x12, v),
  surveys=z12,
  counts=y12)
if (interactive())
  save(josm, file="josm.rda")
