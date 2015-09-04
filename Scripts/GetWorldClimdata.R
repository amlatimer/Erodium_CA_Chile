# Get worldclim / bioclim data for Erodium population locations
# Andrew Latimer 16 April 2012

setwd("/Users/latimer/Documents/Erodium/Greenhouse")

library(maptools); library(fields); library(rgdal); library(gpclib); gpclibPermit(); library(raster); library(gdata)

# use extract() to get values at point locations
#pops = read.xls("popworldclimdata.xlsx")
pops = read.table("additionalpops.csv", sep=",", header=T)

poplocs = SpatialPoints(pops[,c("londd", "latdd")])
#proj4string(poplocs) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

popdata = matrix(0, dim(pops)[1], 19)
for (i in 1:9) {
	wc = raster(paste("/Volumes/LaCie/GISData/Worldclim/bio1-9_30s_bil/bio_", i, ".bil", sep=""))
	popdata[,i] = extract(wc, poplocs)
}

for (i in 10:19) {
	wc = raster(paste("/Volumes/LaCie/GISData/Worldclim/bio10-19_30s_bil/bio_", i, ".bil", sep=""))
	popdata[,i] = extract(wc, poplocs)
}
popdata[popdata>64000] = popdata[popdata>64000]-65536 # correct negative numbers
popdata[,c(1:11)] = popdata[,c(1:11)]/10 # convert to degrees C
popdata = data.frame(popdata)
for (i in 1:19) names(popdata)[i] = paste("BIO", i, sep="")

# Add descriptive columns
popdata = cbind(pops[,c("population", "latdd", "londd")], popdata)

#write.table(popdata, "popworldclimdata.csv", sep=",", row.names=FALSE, col.names=T)
write.table(popdata, "worldclimdata_extra_pops.csv", sep=",", row.names=FALSE, col.names=T)


# check it out

plot(popdata[,4:14])
plot(popdata[,15:22])


chilepops = pops[pops$latdd<0,]

plot(wc, xlim=c(min(chilepops$londd), max(chilepops$londd)), ylim=c(min(chilepops$latdd), max(chilepops$latdd)))
points(poplocs)
plot(pops$latdd[pops$latdd<0], popdata$BIO2[pops$latdd<0])