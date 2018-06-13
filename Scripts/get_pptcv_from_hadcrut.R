# Get global precip and pptcv values from Hadcrut data

setwd("/users/latimer/Google Drive/Erodium_CA_Chile")
library(RNetCDF); library(raster)

popdata = read.csv("./Data/allpopsinfo.csv")
#pops = unique(siteplast.nut$Site)
pops = unique(popdata$Population)
popdata = popdata[match(pops, popdata$Population),]

d = open.nc("./data/cru_ts3.23.2001.2010.pre.dat.nc")

print.nc(d)

# create raster bricks from the netcdf files
r3 = brick("./data/cru_ts3.23.2001.2010.pre.dat.nc")
r2 = brick("./data/cru_ts3.23.1991.2000.pre.dat.nc")
r1 = brick("./data/cru_ts3.23.1981.1990.pre.dat.nc")
nlayers(r1)

# extract monthly precip for each decade, then bind them together
ppt = extract(r1, popdata[,c("Longitude", "Latitude")])
dim(ppt)
sum(!is.na(ppt))
ppt = cbind(ppt, extract(r2, popdata[,c("Longitude", "Latitude")]))
ppt = cbind(ppt, extract(r3, popdata[,c("Longitude", "Latitude")]))
dim(ppt)
sum(!is.na(ppt))
hist(ppt)
hist(log10(ppt[ppt>0])) # no apparent huge outliers

# summarize to annual scale, using "water years" 
# Note water year for CA is Sept-Aug; water year for Chile is March-Feb
ppt = as.data.frame(ppt)
ppt_long = as.data.frame(t(ppt))
names(ppt_long) = pops
ppt_long$year = rep(1981:2010, times=rep(12, 30))
ppt_long$month = rep(1:12, 30)
head(ppt_long)

get_water_year= function(ppt, year, month, wateryears, region) { 
  w = rep(NA, length(wateryears))
  if (region == "CA") {
    for (i in 1:length(wateryears)) {w[i] = sum(ppt[year==wateryears[i] & month<9]) + sum(ppt[year==(wateryears[i]-1) & month>8])}
    return(w)
  }
  if (region=="Chile") {}
    for (i in 1:length(wateryears)) {w[i] = sum(ppt[year==wateryears[i] & month<3]) + sum(ppt[year==(wateryears[i]-1) & month>2])}
   return(w)
}


get_water_year(ppt=ppt_long[,2], year=ppt_long$year, month=ppt_long$month, wateryears=1982:2010, region="CA")
get_water_year(ppt=ppt_long[,1], year=ppt_long$year, month=ppt_long$month, wateryears=1982:2010, region="Chile")

#siteregions = c("Chile", "CA", "Chile", "Chile", "Chile", "Chile", "CA", "Chile", "CA", "Chile", "CA", "CA", "Chile", "CA", "CA", "Chile", "CA", "CA", "Chile")
siteregions = c(rep("CA", 21), rep("Chile",14))
sitecols = sapply(siteregions, FUN=switch, Chile="orange2", CA="royalblue")


ppt_annual=matrix(NA, 29, 35)
for (i in 1:length(siteregions)) {
    ppt_annual[,i] = get_water_year(ppt=ppt_long[,i], year=ppt_long$year, month=ppt_long$month, wateryears=1982:2010, region=siteregions[i])
}
ppt_annual = as.data.frame(ppt_annual)
names(ppt_annual) = pops
head(ppt_annual)

ppt_mean = apply(ppt_annual, 2, mean)
pptcv = apply(ppt_annual, 2, cv)
plot(pptcv~ppt_mean, pch=16, col=sitecols)
# Well, pptcv just IS higher in chile, while rainfall is lower

haddata = data.frame(Site=pops, had_pptmean = ppt_mean, had_pptcv = pptcv)

write.csv(haddata, "hadley_ppt_data_by_site.csv")


# Link the hadley data to greenhouse measurements 

siteplast.water = merge(siteplast.water, haddata, by="Site")
siteplast.nut = merge(siteplast.nut, haddata, by="Site")
d.g.fam = merge(d.g.fam, haddata, by="Site")
sitedata = merge(sitedata, haddata, by="Site")
sitedata_wp = merge(sitedata_wp, haddata, by="Site")
sitedata_np = merge(sitedata_np, haddata, by="Site")
sitedata_lowater = merge(sitedata_lowater, haddata, by="Site")

# Output the data

write.csv(siteplast.water, "siteplast.water.csv")
write.csv(siteplast.nut, "siteplast.nut.csv")
write.csv(d.g.fam, "d.g.fam.csv")
write.csv(sitedata, "sitedata.csv")
write.csv(sitedata_wp, "sitedata_wp.csv")
write.csv(sitedata_np, "sitedata_np.csv")
write.csv(sitedata_lowater, "sitedata_lowater.csv")



