# Explore distributions of traits in field vs greenhouse

library(lattice); library(fields); library(MASS); library(reldist)
setwd("/users/latimer/documents/erodium/alldata2013/")


d = read.csv("All field data combined.csv")

d.g = read.csv("second generation dataset.csv")

z = d$SLA[d$SLA<100000]
par(mfrow=c(2,1))
hist(z/1000, xlim=c(0, 100), nclass=12)
hist(d.g$SLA/1000)

# Functions used:

sitemean <- function(x, site.index) { 
	sites = unique(site.index)
	sites = sites[sites!=""]
	nsites=length(sites)
	x.mean = rep(0, nsites)
	for (i in 1:nsites) x.mean[i] = mean(x[site.index == sites[i]], na.rm=T)
	return(x.mean)
}

sitevar <- function(x, site.index) { 
	sites = unique(site.index)
	sites = sites[sites!=""]
	nsites=length(sites)
	x.mean = rep(0, nsites)
	for (i in 1:nsites) x.mean[i] = var(x[site.index == sites[i]], na.rm=T)
	return(x.mean)
}

sitecv <- function(x, site.index) { 
	return(sitevar(x, site.index)/sitemean(x, site.index))
}


siteplast <- function(trait, site, treatment1, treatment2) {
	return(sitemean(trait[treatment1], site[treatment1]) - sitemean(trait[treatment2], site[treatment2]))
}

plot_plast <- function(reg_water, reg_nut, dataset) {
	plot(reg_water, data=dataset, col = "slateblue",pch=16, cex=2, xlab=as.character(reg_water[3]), ylab=as.character(reg_water[2]))
abline(lm(reg_water, data=dataset)$coef, lwd=2, col="slateblue")
points(reg_nut, data=dataset, col = "orange", pch=16, cex=2)
abline(lm(reg_nut, data=dataset)$coef, lwd=2, col="orange")
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))
}

# make a color scale for treatments
trmtcols = as.character(d.g$treatment)
trmtcols[trmtcols=="high_water"] = "blue"
trmtcols[trmtcols=="low_water"] = "yellow"
trmtcols[trmtcols=="high_nut"] = "green"
trmtcols[trmtcols=="low_nut"] = "orange"
trmtcols.all = as.character(d.g$treatment)
trmtcols.all[trmtcols.all=="high_water"] = "blue"
trmtcols.all[trmtcols.all=="low_water"] = "black"
trmtcols.all[trmtcols.all=="high_nut"] = "green"
trmtcols.all[trmtcols.all=="low_nut"] = "orange"



# Make a data set that summarizes mean traits by family for 
# high water, low water, high nutrients, low nutrients.

ncols = dim(d.g)[2]
numeric.cols = c(1:3, 13:96) # which columns are numeric and thus should be averaged
character.cols = 4:12
for (i in 1:length(numeric.cols)) d.g[,numeric.cols[i]] = as.numeric(d.g[,numeric.cols[i]])
for (i in 1:length(character.cols)) d.g[,character.cols[i]] = as.character(d.g[,character.cols[i]])

families = unique(d.g$field_family_ID)
nfam = length(families)
treatments = unique(d.g$treatment)[1:4]
ntreatments = 4
famdata = as.data.frame(matrix(0, nfam*ntreatments, ncols))
for (i in 1:length(numeric.cols)) famdata[,numeric.cols[i]] = as.numeric(famdata[,numeric.cols[i]])
for (i in 1:length(character.cols)) famdata[,character.cols[i]] = as.character(famdata[,character.cols[i]])
names(famdata) = names(d.g)
for (i in 1:nfam) {
	famindex = d.g$field_family_ID == families[i]
		for (j in 1:ntreatments) {
			trmtindex = d.g$treatment==treatments[j]
			famdata[(i-1)*4+j,numeric.cols] = apply(d.g[famindex & trmtindex, numeric.cols], 2, mean, na.rm=T)	
			z = grep(T, famindex & trmtindex)[1]
			famdata[(i-1)*4+j,character.cols] = d.g[z, character.cols]
	}
}
write.table(famdata, "greenhouse_data_by_family.csv", sep=",", row.names=F, col.names=T)



##### Explore at what levels of grouping the trait means are varying (among family within pop, among pop, among region)

famdata = read.csv("greenhouse_data_by_family.csv")

m = lm(log(days_to_flower)~region+Site, famdata, subset=famdata$treatment=="low_nut")
summary(m); anova(m)
# Most of the variance (~80%) is explained by source population. Strong local genetic differences in this critical trait. 

boxplot(log(days_to_flower)~Site, famdata[famdata$treatment=="low_nut",], notch=T)

# Is it connected to fitness?
summary(lmer(log(gen2_total_seeds+0.1)~days_to_flower + (1+days_to_flower|Site), famdata[famdata$treatment=="low_nut",]))
# yes -- and apparently no within-site variation in slope, only intercept.
xyplot(gen2_total_seeds~days_to_flower|Site, famdata[famdata$treatment=="low_nut",])
plot(gen2_total_seeds~days_to_flower, famdata[famdata$treatment=="low_nut", ])
text(famdata$days_to_flower[famdata$treatment=="low_nut"], famdata$gen2_total_seeds[famdata$treatment=="low_nut"], label = famdata$Site [famdata$treatment=="low_nut"], cex=0.5)
# Yes, but most of the selection is among the populations, not among families within populations. 

###### Now check for plasticity of the traits
famdata = famdata[!is.na(famdata$field_family_ID),]
famplastdata = famdata[famdata$treatment %in% c("high_water"),]

numeric.cols = c(1:3, 13:97) # which columns are numeric and thus should be differenced

families = unique(d.g$field_family_ID)

# remove the families not represented in one of the treatments so can't have plasticity value. 
families = families[-c(53, 68, 96, 130, 133, 152, 166, 168)]
famplastdata = famplastdata[famplastdata$field_family_ID %in% families,]


nfam = length(families)
treatments = unique(d.g$treatment)[1:4]
waterdiff_days_to_flower = rep(NA, nfam)
nutdiff_days_to_flower = rep(NA, nfam)
waterdiff_length_longest_stem = rep(NA, nfam)
nutdiff_length_longest_stem = rep(NA, nfam)
waterdiff_mean_seeds = rep(NA, nfam)
nutdiff_mean_seeds = rep(NA, nfam)

for (i in 1:nfam) {
	famindex = famdata$field_family_ID == families[i]
	waterdiff_length_longest_stem[i] = famdata$gen2_length_longest_stem[famindex & famdata$treatment=="high_water"] - famdata$gen2_length_longest_stem[famindex & famdata$treatment=="low_water"]
	nutdiff_length_longest_stem[i] = famdata$gen2_length_longest_stem[famindex & famdata$treatment=="high_nut"] - famdata$gen2_length_longest_stem[famindex & famdata$treatment=="low_nut"]
	waterdiff_days_to_flower[i] = famdata$days_to_flower[famindex & famdata$treatment=="high_water"] - famdata$days_to_flower[famindex & famdata$treatment=="low_water"]
	nutdiff_days_to_flower[i] = famdata$days_to_flower[famindex & famdata$treatment=="high_nut"] - famdata$days_to_flower[famindex & famdata$treatment=="low_nut"]
	waterdiff_mean_seeds[i] = mean(famdata$days_to_flower[famindex & famdata$treatment %in% c("high_nut", "low_nut")], na.rm=T)
	nutdiff_mean_seeds[i] = mean(famdata$days_to_flower[famindex & famdata$treatment %in% c("high_water", "low_water")], na.rm=T)
}
	famplastdata = cbind(famplastdata, waterdiff_length_longest_stem, waterdiff_days_to_flower, waterdiff_mean_seeds, nutdiff_length_longest_stem, nutdiff_days_to_flower, nutdiff_mean_seeds)

# Add variance and CV in site nutrients and cover
d$population = as.character(d$population)
sites = unique(d.g$Site)
nsites = length(sites)
var_total_cover = rep(NA, nsites)
mean_total_cover = rep(NA, nsites)
var_P1 = rep(NA, nsites)
mean_P1 = rep(NA, nsites)
var_soil_water = rep(NA, nsites)
mean_soil_water = rep(NA, nsites)

for (i in 1:nsites) {
	z = d$population == sites[i]
	var_total_cover[i] = var(d$total[z], na.rm=T)
	mean_total_cover[i] = mean(d$total[z], na.rm=T)
	var_P1[i] = var(d$P1[z], na.rm=T)
	mean_P1[i] = mean(d$P1[z], na.rm=T)
	var_soil_water[i] = var(d$soil_water_content[z], na.rm=T)
	mean_soil_water[i] = mean(d$soil_water_content[z], na.rm=T)
}
z = match(famplastdata$Site, sites)
famplastdata$var_total_cover = var_total_cover[z]
famplastdata$mean_total_cover = mean_total_cover[z]
famplastdata$var_P1 = var_P1[z]
famplastdata$mean_P1 = mean_P1[z]
famplastdata$var_soil_water = var_soil_water[z]
famplastdata$mean_soil_water = mean_soil_water[z]

write.table(famplastdata, "greenhouse_plasticity_data_by_family.csv", sep=",", row.names=F, col.names=T)
	
	
famplastdata=read.table("greenhouse_plasticity_data_by_family.csv", sep=",", header=T)

	
# PLots vs source site precipitation
par(mfrow=c(1,2))
plot(waterdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="Precip", main="California", ylab="Plasticity in flowering")
abline(lm(waterdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="California",])$coef, lwd=2, col="slateblue")
points(nutdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="California",], col = "orange", pch=16, cex=2)
abline(lm(nutdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="California",])$coef, lwd=2, col="orange")

plot(waterdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="Precip", main="Chile", ylab="Plasticity in flowering")
abline(lm(waterdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="Chile",])$coef, lwd=2, col="slateblue")
points(nutdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
abline(lm(nutdiff_days_to_flower~BIO12, famplastdata[famplastdata$region=="Chile",])$coef, lwd=2, col="orange")
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))

# Very similar negative relationship in both regions between plasticity in first flower date and source site precipitation. higher precip = greater maladaptive plasticity in flower date. 


# Plot vs interannual and intra-annual variation in precip 
famplastdata$drywetratio = famplastdata$BIO17/famplastdata$BIO16

plot(waterdiff_days_to_flower~drywetratio, famplastdata, col = "slateblue",pch=16, cex=2, xlab="Precip", ylab="Plasticity in flowering")
abline(lm(waterdiff_days_to_flower~precipCV, famplastdata[famplastdata$region=="California",])$coef, lwd=2, col="slateblue")
points(nutdiff_days_to_flower~precipCV, famplastdata[famplastdata$region=="California",], col = "orange", pch=16, cex=2)
abline(lm(nutdiff_days_to_flower~precipCV, famplastdata[famplastdata$region=="California",])$coef, lwd=2, col="orange")

plot(precipCV~BIO12, famplastdata)
plot((BIO17/BIO16)~BIO12, famplastdata)
# Nothing much here????


# PLots vs source site variance in soil nutrients (phosphorus)
# Nothing dramatic here as far as I can tell. 




plot_plast(waterdiff_days_to_flower~var_total_cover, nutdiff_days_to_flower~var_total_cover, dataset=famplastdata)

plot_plast(waterdiff_days_to_flower~var_P1, nutdiff_days_to_flower~var_P1, dataset=famplastdata)
plot_plast(waterdiff_length_longest_stem~var_P1, nutdiff_length_longest_stem~var_P1, dataset=famplastdata)
plot_plast(waterdiff_length_longest_stem~mean_total_cover, nutdiff_length_longest_stem~mean_total_cover, dataset=famplastdata)
plot_plast(waterdiff_length_longest_stem~var_total_cover, nutdiff_length_longest_stem~var_total_cover, dataset=famplastdata)

# compare regions
par(mfrow=c(1,2))
plot_plast(waterdiff_days_to_flower~var_P1, nutdiff_days_to_flower~var_P1, dataset=famplastdata[famplastdata$region=="California",])
plot_plast(waterdiff_days_to_flower~var_P1, nutdiff_days_to_flower~var_P1, dataset=famplastdata[famplastdata$region=="Chile",])

#compare means vs plasticity effects
par(mfrow=c(1,2))
plot_plast(waterdiff_days_to_flower~mean_P1, nutdiff_days_to_flower~var_P1, dataset=famplastdata)
plot_plast(waterdiff_days_to_flower~var_P1, nutdiff_days_to_flower~var_P1, dataset=famplastdata)

# Really not a lot here. Maybe plasticity in days to flower in response to nutrients in California only is related to site heterogeneity -- I think this is Brooke's finding. 


# Compare field to greenhouse values for SLA -- Nothing there. 
families = unique(d.g$field_family_ID)
z = match(families, d$ID)
d.g.temp = d.g[d.g$treatment=="low_nut",]
z.g = match(families, d.g.temp$field_family_ID)
plot(d$SLA[z], d.g.temp$SLA[z.g], xlim=c(0, 50000))

plot(log(SLA)~days_to_flower, d.g, col=factor(d.g$region))#trmtcols)

# DO the same at the population level. 

sites = unique(d.g$Sites)
nsites - length(sites)



### Exploration for regression class:
d = read.table("greenhouse_data_by_family.csv", sep=",", header=T)

d = d[d$treatment %in% c("high_water", "low_water"),]
d$gen2_length_longest_stem[d$gen2_length_longest_stem==0] = NA


summary(lm(log(gen2_length_longest_stem)~days_to_flower*treatment, d))

plot(log(gen2_length_longest_stem)~days_to_flower, d, pch=16, col=(d$treatment=="low_water")+1)

# Q how would you interpret the intercept and the effect of treatment? (remember the interaction term)
# Does this give a sensible interpretation? 
      
# One reason to center the explanatory variables is to make the effects more interpretable. 
      
d$days_to_flower=d$days_to_flower-mean(d$days_to_flower, na.rm=T)      
      
summary(lm(log(gen2_length_longest_stem)~days_to_flower*treatment, d))
      
plot(log(gen2_length_longest_stem)~days_to_flower, d, pch=16, col=(d$treatment=="low_water")+1)

# Now how to interpret the intercept and treament effects? 
    
# Q how would you compare the sizes of the effect of treatment and days to flower? 
# G&H recommend centering all variables, and dividing continuous variables by 2x their standard deviation. 
# This makes the effects of continuous variables more comparable to those of 1/0 categorical variables. 
# See G&H pp56-57.

d.scaled = d

plot(days_to_flower~gen2_length_longest_stem, d, pch=16, col=(d$region=="California")+1)

summary(lm(days_to_flower~gen2_length_longest_stem*treatment, d)# subset=d$region=="Chile"))

summary(lm(days_to_flower~SLA*treatment, d)#, subset=d$region=="Chile"))

summary(lm(days_to_flower~BIO12*treatment, d, subset=d$region=="Chile"))
summary(lm(days_to_flower~BIO12*region, d))

xyplot((days_to_flower+90)~BIO12|treatment, d, xlab="Mean annual precip (mm)", ylab="Days to flower", cex.axis=1.2, cex.lab=1.5)

# Aggregate to pop level
d_pop = aggregate(d[,c("days_to_flower", "BIO12")], by=list(d$Site), FUN=mean)
head(d_pop)
plot((days_to_flower+90)~BIO12, d_pop, xlab="Mean annual precip (mm)", ylab="Days to flower", cex.axis=1.2, cex.lab=1.5)
  
  
  
palette(tim.colors(round(max(famplastdata$SLA, na.rm=T)/1000)))
xyplot(waterdiff_days_to_flower~BIO12|treatment, famplastdata, col=famplastdata$SLA/1000, pch=16)


palette(tim.colors(round(max(famplastdata$days_to_flower, na.rm=T))))
xyplot(waterdiff_days_to_flower~SLA, famplastdata, col=famplastdata$days_to_flower, pch=16)

  
xyplot(waterdiff_days_to_flower~SLA|treatment, famplastdata)#, col=famplastdata$SLA/1000)

xyplot(waterdiff_days_to_flower~days_to_flower, famplastdata, col=famplastdata$SLA/1000)#, col=famplastdata$SLA/1000)
  
  
xyplot(log(SLA)~BIO12|treatment, famplastdata)#, col=famplastdata$SLA/1000)
  
summary(lm(waterdiff_days_to_flower~BIO12+log(SLA), famplastdata))
  
# Check correlation between source site precip and cover
pairs(d[,c("BIO12", "grass_", "forb", "total_cover")])
# nothing there!
