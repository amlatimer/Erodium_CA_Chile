# Explore distributions of reproductive effort (fitness proxy) in field populations
library(lattice); library(fields); library(MASS); library(reldist)
setwd("/users/latimer/documents/erodium/alldata2013/")


d = read.csv("All field data combined.csv")

d.g = read.csv("second generation dataset.csv")

d.g$length_per_internode = d.g$gen2_length_longest_stem/d.g$gen2_no_internodes

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


# Get index of populations

pops = unique(d$population)
npops = length(pops)

# only pops used in experiment
ghousepops = as.character(unique(d.g$Site)[1:20])

# Q how many indivs per pop? 
table(d$population)


# Calculate Gini coefs for these pops
gini.pop = rep(0, 20)
cv.pop = rep(0, 20)
negbin.pop = rep(0, 20)
for (i in 1:20) {
	y = d$reproductive_effort[d$population==ghousepops[i]]
	gini.pop[i] = gini(y)
	cv.pop[i] = var(y)/mean(y)
	negbin.pop[i] = as.numeric(fitdistr(y, "negative binomial")[[1]][1])
}

fitness.gini = data.frame(ghousepops, gini.pop)
write.table(fitness.gini, "fitness.gini.csv", sep=",", row.names=F, col.names=T)

# plot distribution of fitness for each pop
xyplot(log(reproductive_effort)~grass|population, data=d)

par(mfrow=c(5, 4), mar=rep(1.1, 4))
for (i in 1:20) {
	y = d$reproductive_effort[d$population==ghousepops[i]]
	qqplot(x=rpois(1000, lambda=mean(y)), y, plot.it=T, main = round(gini.pop[i],2))
}




# Explore greenhouse results for elongation 

xyplot(gen2_length_longest_stem~total_cover|Site, d.g)
xyplot(log(length_per_internode)~total_cover|Site, d.g)


calindex = c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)

trait = sitemean(d.g$gen2_length_longest_stem, d.g$Site)
envvar = sitemean(d.g$grass_, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1)


plot(d.g[,39:42])
cor(d.g[,39:42], use="pairwise.complete")

trait = sitemean(d.g$gen2_length_longest_stem, d.g$Site)
envvar = sitemean(d.g$Erodium_total, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1)
text(envvar, trait, sites, trait,col = calindex+1)

plot(forb~grass_, d.g)
plot(forb~grass_, d.g)

par(mfrow=c(3, 2))
plot(log(gen2_length_longest_stem)~grass_, d.g[d.g$region=="California",], main="California")
plot(log(gen2_length_longest_stem)~grass_, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(log(gen2_length_longest_stem)~forb, d.g[d.g$region=="California",], main="California")
plot(log(gen2_length_longest_stem)~forb, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(log(gen2_length_longest_stem)~Erodium_total, d.g[d.g$region=="California",], main="California")
plot(log(gen2_length_longest_stem)~Erodium_total, d.g[d.g$region=="Chile",], col=2, main="Chile")


par(mfrow=c(3, 2))
plot(gen2_no_internodes~grass_, d.g[d.g$region=="California",], main="California")
plot(gen2_no_internodes~grass_, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(gen2_no_internodes~forb, d.g[d.g$region=="California",], main="California")
plot(gen2_no_internodes~forb, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(gen2_no_internodes~Erodium_total, d.g[d.g$region=="California",], main="California")
plot(gen2_no_internodes~Erodium_total, d.g[d.g$region=="Chile",], col=2, main="Chile")

# Elongation -- calculate length per internode

plot(length_per_internode~gen2_length_longest_stem, d.g)	

par(mfrow=c(4, 2))
plot(log(length_per_internode)~total_cover, d.g[d.g$region=="California",], main="California")
plot(log(length_per_internode)~total_cover, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(log(length_per_internode)~grass_, d.g[d.g$region=="California",], main="California")
plot(log(length_per_internode)~grass_, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(log(length_per_internode)~forb, d.g[d.g$region=="California",], main="California")
plot(log(length_per_internode)~forb, d.g[d.g$region=="Chile",], col=2, main="Chile")
plot(log(length_per_internode)~Erodium_total, d.g[d.g$region=="California",], main="California")
plot(log(length_per_internode)~Erodium_total, d.g[d.g$region=="Chile",], col=2, main="Chile")

summary(glm(log(length_per_internode)~Erodium_total, data=d.g, subset=d.g$region=="Chile"))
summary(lm(length_per_internode~Erodium_total, data=d.g, subset=d.g$region=="California"))

summary(glm(log(length_per_internode)~grass_, data=d.g, subset=d.g$region=="Chile"))
summary(glm(log(length_per_internode)~grass_, data=d.g, subset=d.g$region=="California"))

summary(glm(log(length_per_internode)~forb, data=d.g, subset=d.g$region=="Chile"))
summary(glm(log(length_per_internode)~forb, data=d.g, subset=d.g$region=="California"))


summary(glm(log(length_per_internode)~Erodium_total+grass_+forb, data=d.g, subset=d.g$region=="Chile"))
summary(glm(log(length_per_internode)~Erodium_total*grass_*forb, data=d.g, subset=d.g$region=="California"))

summary(glm(log(gen2_length_longest_stem+1)~Erodium_total+grass_+forb, data=d.g, subset=d.g$region=="Chile"))
summary(glm(log(gen2_length_longest_stem+1)~Erodium_total*grass_*forb, data=d.g, subset=d.g$region=="California"))

lmer(log(gen2_length_longest_stem+1)~Erodium_total+grass_+forb, data=d.g, subset=d.g$region=="California")

plot(d.g[d.g$region=="California",39:41])
plot(d.g[d.g$region=="Chile",39:41])

# some plots at the population level

trait = sitemean(d.g$gen2_length_longest_stem, d.g$Site)
envvar = sitemean(d.g$Erodium_total, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1)
text(envvar, trait, sites, trait,col = calindex+1)

summary(lm(log(trait[1:10])~envvar[1:10]))
summary(lm(log(trait[11:20])~envvar[11:20]))

trait = sitemean(d.g$length_per_internode, d.g$Site)
envvar = sitemean(d.g$Erodium_total, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1, xlab="Erodium total cover", ylab = "Stem length per internode", cex.lab=1.5, cex.axis=1.25)
text(envvar, trait+0.05, sites, trait,col = calindex+1)

trait = sitemean(d.g$length_per_internode, d.g$Site)
envvar = sitemean(d.g$BIO12, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1)
text(envvar, trait, sites, trait,col = calindex+1)

# Q what environmental variables are associated with Erodium cover in CA?

plot(Erodium_total~management.type, d.g[d.g$region=="California",])
plot(Erodium_total~P1, d.g[d.g$region=="California",]) # inverse rel of P to Erodium cover. 
plot(log(length_per_internode)~P1, d.g[d.g$region=="California",]) 
# inverse rel of P to Erodium cover, but P isn't associated with elongation. 

plot(grass_~P1, d.g[d.g$region=="California",]) # surprisingly, negative assn.
plot(Erodium_total~grass_, d.g[d.g$region=="California",])
plot(length_per_internode~grass_, d.g[d.g$region=="California",]) 


# show relationship of internode length vs stem length
plot(length_per_internode~gen2_length_longest_stem, d.g[d.g$region=="California",], col = d.g$Site)
# check for variation among sites in this relationship
elong.slopes = as.numeric(coef(lm(length_per_internode~gen2_length_longest_stem %in% Site, d.g, subset = d.g$region=="California"))[2:11])
envvar = sitemean(d.g$BIO12, d.g$Site)
plot(elong.slopes~envvar[1:10])
# Doesn't seem to be much here!

boxplot(Erodium_total~climate, d.g[d.g$region=="California",])
boxplot(log(length_per_internode)~climate, d.g[d.g$region=="California",])
# Intriguingly, Erodium cover generally higher inland, but this does not seem to be what's driving association of elongation with Erodium cover. Little association of elongation with inland/coastal location. 

# BUT: Precipitation at a site explains a lot of site-level elongation
trait = sitemean(d.g$length_per_internode, d.g$Site)
envvar = sitemean(d.g$BIO12, d.g$Site)
plot(trait~envvar, pch=16, cex=1.5, col = calindex+1)
text(envvar, trait, sites, trait,col = calindex+1)

plot(log(length_per_internode)~BIO12, d.g)
plot(log(gen2_length_longest_stem)~BIO12, d.g)
plot(Erodium_total~BIO12, d.g) # weirdly Erodium total cover not super strongly associated with precip, even though source-site Precip is. 

# How much variation in stem length can source site precip and Erodium cover together explain? 
summary(lm(log(length_per_internode)~BIO12*Erodium_total, d.g))
# At the population level:
#trait = sitemean(d.g$length_per_internode, d.g$Site)
#trait = sitemean(d.g$gen2_length_longest_stem, d.g$Site)
trait = sitemean(d.g$gen2_no_internodes, d.g$Site)
envvar1 = sitemean(d.g$BIO12, d.g$Site)
envvar2 = sitemean(d.g$Erodium_total, d.g$Site)
summary(lm(log(trait)~envvar1*envvar2))

# SO: plant stem length and number of internodes is strongly associated with source-site precipitation, but not with Erodium cover
plot(sitemean(d.g$BIO12, d.g$Site), sitemean(d.g$gen2_no_internodes, d.g$Site), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = "Number of internodes in greenhouse")
legend("topright",c("Chile", "California"), col=c(1, 2), pch=c(16, 16))
# Whereas "elongation" in the sense of the length of internodes is strongly associated with source-site Erodium cover, but not with precipitation. 
plot(sitemean(d.g$Erodium_total, d.g$Site), sitemean(d.g$length_per_internode, d.g$Site), col = calindex+1, pch=16, xlab="Erodium cover", ylab = "Stem length per internode")
legend("bottomright",c("Chile", "California"), col=c(1, 2), pch=c(16, 16))

# Check whether these results are contingent on greenhouse treatment. 
par(mfrow=c(2,2))
trmts = unique(d.g$treatment)
for (i in 1:4) {
	dtemp = d.g[d.g$treatment==trmts[i],]
	plot(sitemean(dtemp$BIO12, dtemp$Site), sitemean(dtemp$gen2_no_internodes, dtemp$Site), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Number of internodes in greenhouse", main=trmts[i])
}
# Yes this one is remarkably consistent, even with indiv pops being obvious across treatments. Strong genetic effect. 
par(mfrow=c(2,2))
trmts = unique(d.g$treatment)
for (i in 1:4) {
	dtemp = d.g[d.g$treatment==trmts[i],]
	plot(sitemean(dtemp$Erodium_total, dtemp$Site), sitemean(dtemp$length_per_internode, dtemp$Site), col = calindex+1, pch=16, xlab="Erodium cover", ylab = "Stem length per internode", main=trmts[i])
	#legend("bottomright",c("Chile", "California"), col=c(1, 2), pch=c(16, 16))}
}
# Yes these too, though in the low-water treamtent, Kearney jumps out as an outlier that breaks the pattern. Generally no pattern in the Chilean sites, only CA. 


# Now check plasticity in these traits with respect to water treatment.
siteplast <- function(trait, site, treatment1, treatment2) {
	return(sitemean(trait[treatment1], site[treatment1]) - sitemean(trait[treatment2], site[treatment2]))
}

par(mfrow=c(2,2))
plot(sitemean(d.g$BIO12, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_water", d.g$treatment=="low_water"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in # internodes", main="High vs low water")

plot(sitemean(d.g$BIO12, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_water", d.g$treatment=="low_water"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in internode length", main="High vs low water")

plot(sitemean(d.g$Erodium_total, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_water", d.g$treatment=="low_water"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in # internodes", main="High vs low water")

plot(sitemean(d.g$Erodium_total, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_water", d.g$treatment=="low_water"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in internode length", main="High vs low water")

# Not much there, except in internode length vs water treatment in CA only. 

# Same but now check plasticity with respect to nutrient treatement
par(mfrow=c(2,2))
plot(sitemean(d.g$OM, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in # internodes", main="High vs low nutrients")

plot(sitemean(d.g$OM, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in internode length", main="High vs low nutrients")

plot(sitemean(d.g$OM, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in # internodes", main="High vs low nutrients")

plot(sitemean(d.g$OM, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in internode length", main="High vs low nutrients")

# Yes higher plasticity in internode length is (noisily) associated with drier source sites in both regions -- so stems are longer on average, and their internode lengths are also more plastic, in the dry populations. 

# Q whether plasticity in stem length might be related to local heterogeneity

par(mfrow=c(2,2))
plot(sitecv(d.g$P1, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in # internodes", main="High vs low nutrients")

plot(sitecv(d.g$P1, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Mean Annual Precip (mm)", ylab = 	"Plasticity in internode length", main="High vs low nutrients")

plot(sitecv(d.g$P1, d.g$Site), siteplast(d.g$gen2_no_internodes, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in # internodes", main="High vs low nutrients")

plot(sitecv(d.g$P1, d.g$Site), siteplast(d.g$length_per_internode, d.g$Site, d.g$treatment=="high_nut", d.g$treatment=="low_nut"), col = calindex+1, pch=16, xlab="Source site Erodium cover", ylab = 	"Plasticity in internode length", main="High vs low nutrients")





# Check if source site precip is correlated with nutrients and Erodium cover
plot(d.g[,c("BIO12", "P1", "HCO3_P", "OM", "Erodium_total")])
cor(d.g[,c("BIO12", "P1", "HCO3_P", "OM", "Erodium_total")], use="pairwise.complete")
# Associated moderately with OM, weakly with P

# Standardize variables for this analysis
standardize <- function(x) { return((x-mean(x, na.rm=T))/sd(x, na.rm=T)) }
d.g[,c("BIO12", "P1", "HCO3_P", "OM", "Erodium_total")] = apply(d.g[,c("BIO12", "P1", "HCO3_P", "OM", "Erodium_total")], 2, standardize)

summary(lm(length_per_internode~Erodium_total+BIO12+P1+OM, d.g)) # The best association by far is with Erodium cover, then precip
summary(lm(gen2_no_internodes~Erodium_total+BIO12+P1+OM, d.g)) # Here all 4 factors are associated: Erodium cover, precip, P and OM.

# NEXT: Need to do this analysis for the family means and plasticities. 

# Make a data set that summarizes greenhouse data to the population level. 



ncols = dim(d.g)[2]
numeric.cols = c(1:3, 13:97) # which columns are numeric and thus should be averaged
character.cols = 4:12
for (i in 1:length(numeric.cols)) d.g[,numeric.cols[i]] = as.numeric(d.g[,numeric.cols[i]])
for (i in 1:length(character.cols)) d.g[,character.cols[i]] = as.character(d.g[,character.cols[i]])

groups = unique(d.g$Site)
ngroups = length(groups)
treatments = unique(d.g$treatment)[1:4]
ntreatments = 4
groupdata = as.data.frame(matrix(0, ngroups*ntreatments, ncols))
for (i in 1:length(numeric.cols)) groupdata[,numeric.cols[i]] = as.numeric(groupdata[,numeric.cols[i]])
for (i in 1:length(character.cols)) groupdata[,character.cols[i]] = as.character(groupdata[,character.cols[i]])
names(groupdata) = names(d.g)
for (i in 1:ngroups) {
	groupindex = d.g$Site == groups[i]
		for (j in 1:ntreatments) {
			trmtindex = d.g$treatment==treatments[j]
			groupdata[(i-1)*4+j,numeric.cols] = apply(d.g[groupindex & trmtindex, numeric.cols], 2, mean, na.rm=T)	
			z = grep(T, groupindex & trmtindex)[1]
			groupdata[(i-1)*4+j,character.cols] = d.g[z, character.cols]
	}
}

#write.table(groupdata, "greenhouse_by_siteandtreatment.csv", sep=",", row.names=F, col.names=T)


######## ANALYZE POPULATION X TREATMENT LEVEL DATA
groupdata = read.table("greenhouse_by_siteandtreatment.csv", sep=",", header=T)


summary(lm(length_per_internode~gen2_no_internodes+Erodium_total+BIO12+P1+OM, groupdata)) # The best association by far is with Erodium cover, then precip
summary(lm(gen2_no_internodes~Erodium_total+BIO12+P1+OM, groupdata)) # Here all 4 factors are associated: Erodium cover, precip, P and OM.

plot(groupdata[,c("length_per_internode", "Erodium_total", "BIO12", "P1", "OM")], col = as.factor(groupdata$treatment))
# make a color scale for treatments
trmtcols = as.character(groupdata$treatment)
trmtcols[trmtcols=="high_water"] = "blue"
trmtcols[trmtcols=="low_water"] = "yellow"
trmtcols[trmtcols=="high_nut"] = "green"
trmtcols[trmtcols=="low_nut"] = "orange"
trmtcols.all = as.character(d.g$treatment)
trmtcols.all[trmtcols.all=="high_water"] = "blue"
trmtcols.all[trmtcols.all=="low_water"] = "black"
trmtcols.all[trmtcols.all=="high_nut"] = "green"
trmtcols.all[trmtcols.all=="low_nut"] = "orange"


par(mfrow=c(2,2))
plot(length_per_internode~Erodium_total, groupdata, col = trmtcols,pch=16, cex=1.5)
plot(length_per_internode~BIO12, groupdata, col = trmtcols,pch=16, cex=1.5)
points(length_per_internode~BIO12, groupdata[groupdata$region=="California",], col = 1, pch=1, cex=1.5)
plot(length_per_internode~P1, groupdata, col = trmtcols,pch=16, cex=1.5)
plot(length_per_internode~gen2_no_internodes, groupdata, col = trmtcols, pch=16, cex=1.5)
points(length_per_internode~gen2_no_internodes, groupdata[groupdata$region=="California",], col = 1, pch=1, cex=1.5)


par(mfrow=c(2,2))
plot(gen2_no_internodes~Erodium_total, groupdata, col = trmtcols,pch=16, cex=1.5)
plot(gen2_no_internodes~BIO12, groupdata, col = trmtcols,pch=16, cex=1.5)
plot(gen2_no_internodes~P1, groupdata, col = trmtcols,pch=16, cex=1.5)
plot(gen2_no_internodes~length_per_internode, groupdata, col = trmtcols, pch=16, cex=1.5)
points(gen2_no_internodes~length_per_internode, groupdata[groupdata$region=="California",], col = 1, pch=1, cex=1.5)


# Separate the regions to check on relationships with important variables
par(mfrow=c(1,2))
plot(length_per_internode~BIO12, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="California")
plot(length_per_internode~BIO12, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="Chile")
legend("topright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))

par(mfrow=c(1,2))
plot(length_per_internode~Erodium_total, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Erodium cover", main="California")
plot(length_per_internode~Erodium_total, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Erodium cover", main="Chile")
legend("topright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))

par(mfrow=c(1,2))
plot(length_per_internode~P1, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Phosphorus -- P1", main="California")
plot(length_per_internode~P1, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Phosphorus -- P1", main="Chile")
legend("topleft", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))

par(mfrow=c(1,2))
plot(gen2_no_internodes~BIO12, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="California")
plot(gen2_no_internodes~BIO12, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="Chile")
legend("topright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))

par(mfrow=c(1,2))
plot(gen2_no_internodes~Erodium_total, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Erodium cover", main="California")
legend("topleft", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))
plot(gen2_no_internodes~Erodium_total, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Erodium cover", main="Chile")

par(mfrow=c(1,2))
plot(gen2_no_internodes~P1, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Phosphorus -- P1", main="California")
plot(gen2_no_internodes~P1, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Phosphorus -- P1", main="Chile")
legend("bottomright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))






par(mfrow=c(1,2))
plot(days_to_flower~BIO12, groupdata[groupdata$region=="California",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="California")
plot(days_to_flower~BIO12, groupdata[groupdata$region=="Chile",], col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="Chile")
legend("topright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))


# What these show: 
# 1) Stem length (number of internodes) responds consistently to both treatments, and these responses are consistent among populations and regions: In high water and nigh nutrient treatments, long stems, in low treatments, short stems. 

# 2) Sites from high-precipitation sites, and low-Erodium cover sites, have shorter stems across all treatments. There is no strong trend in plasticity along among-population gradients in environment. Outlier: Sweeney Corral -- lowest precip, low elongation and number of internodes. 


# Q1 -- is the higher growth that's associated with source-site precip a function of growing-season temperature, growing season length (and/or flowering time)? 

# A: Worldclim temperature variables not as correlated with response as precip. 
summary(lm(gen2_no_internodes~BIO18*treatment, groupdata))
# BIO16, BIO13, BIO12 are the strongly significant ones, also highly correlated to the point of being indistinguishable (0.99). 
cor(groupdata[,c("BIO12", "BIO13", "BIO16")])
cor(groupdata[,c("gen2_no_internodes", "length_per_internode", "BIO12", "Erodium_total", "P1")])

# Q: to what extent are the patterns in length_per_internode simply driven by its correlation with number of internodes? 
summary(lm(length_per_internode~gen2_no_internodes+Erodium_total, groupdata))
summary(lm(length_per_internode~gen2_no_internodes+BIO12+Erodium_total, groupdata[groupdata$region=="California",]))
summary(lm(length_per_internode~gen2_no_internodes+BIO12+Erodium_total, groupdata[groupdata$region=="Chile",]))
summary(lm(length_per_internode~gen2_no_internodes+BIO12+Erodium_total, groupdata[groupdata$treatment %in% c("high_water", "high_nut"),]))
summary(lm(length_per_internode~gen2_no_internodes+BIO12+Erodium_total, groupdata[groupdata$treatment %in% c("low_water", "low_nut"),]))
# A: Entirely! So length_per_internode doesn't seem to be an independent variable at all... Too bad. 

# Growing season length is presumably associated with precip, which will tend to drive longer growing time. But in greenhouse the plants from dry sites grow much longer stems. So inverse relationship with growing season length in source site.  
boxplot(days_to_flower~treatment, data=groupdata, notch=T)
summary(lm(days_to_flower~gen2_no_internodes+BIO12+Erodium_total+P1, groupdata))
summary(lm(gen2_no_internodes~days_to_flower+BIO12+Erodium_total+P1, groupdata))
cor(groupdata[,c("gen2_ffdate", "gen2_no_internodes")])

# Intriguingly, number of internodes is strongly negatively correlated with first flower date. Taking this correlated growth pattern into account, the effect of precipitation and Erodium cover keep the same sign, while their effects are weakened. So this "syndrome" = rapid growth, long stem, early flowering, is mostly a single reponse, and it's related, somewhat weakly, to source site precip and Erodium cover. 


plot(days_to_flower~gen2_no_internodes, groupdata, col=trmtcols, pch=16, cex=2)

# Make a better color scheme for showing precip
rain.colors = colorRampPalette(c("yellow","cyan","blue"))

# Explore relationship of precip to stem length and flowering time

par(mfrow=c(2,2))
palette(rain.colors(1000))
plot(days_to_flower~gen2_length_longest_stem, groupdata[groupdata$treatment=="high_water",], col=groupdata$BIO12[groupdata$treatment=="high_water"], pch=16, cex=1.5, main="High water", ylim=c(0, 150), xlim=c(1, 20), xlab="Length of longest stem", ylab="First flower date")
plot(days_to_flower~gen2_length_longest_stem, groupdata[groupdata$treatment=="low_water",], col=groupdata$BIO12[groupdata$treatment=="low_water"], pch=16, cex=1.5, main="Low water", ylim=c(0, 150), xlim=c(1, 20), xlab="Length of longest stem", ylab="First flower date")
plot(days_to_flower~gen2_length_longest_stem, groupdata[groupdata$treatment=="high_nut",], col=groupdata$BIO12[groupdata$treatment=="high_nut"], pch=16, cex=1.5, main="High nutrients", ylim=c(0, 150), xlim=c(1, 20), xlab="Length of longest stem", ylab="First flower date")
plot(days_to_flower~gen2_length_longest_stem, groupdata[groupdata$treatment=="low_nut",], col=groupdata$BIO12[groupdata$treatment=="low_nut"], pch=16, cex=1.5, main="Low nutrients", ylim=c(0, 150), xlim=c(1, 20), xlab="Length of longest stem", ylab="First flower date")
legend("topright", c("10", "250", "500", "1000"), col=c(10, 250,500,1000), pch=15, title="Precip (mm)")

palette(rain.colors(max(groupdata$gen2_no_seeds)))
cloud(days_to_flower~BIO12*gen2_length_longest_stem, data=groupdata, colorkey=T, col = groupdata$gen2_no_seeds)



cloud(days_to_flower~BIO12*gen2_length_longest_stem, data=groupdata[z,], colorkey=T, col = groupdata$days_to_flower[z],screen = list(z = 30, x = -30, y = 0))

# Days to flower versus source-site precip, by treatment with fitness color scheme
par(mfrow=c(2,2))
palette(rain.colors(250)); 
z = groupdata$treatment=="high_water"
plot(days_to_flower~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5)
z = groupdata$treatment=="low_water"
plot(days_to_flower~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5)
z = groupdata$treatment=="high_nut"
plot(days_to_flower~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5)
z = groupdata$treatment=="low_nut"
plot(days_to_flower~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5)
legend("bottomright", c("10", "50", "150", "250"), col=c(10, 50, 150,250), pch=15, title="Fitness (seeds)")

# Most populatoins flower between 90 and 100 days, a bit longer in the low-resource treatments. A few populations flower substantially earlier, and these are all from dry sites and all have high fitness.

# Stem length versus source-site precip, by treatment with fitness color scheme
par(mfrow=c(2,2))
palette(rain.colors(250)); 
z = groupdata$treatment=="high_water"
plot(gen2_length_longest_stem~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5, ylim = c(0, 40), main="High water")
z = groupdata$treatment=="low_water"
plot(gen2_length_longest_stem~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5, ylim = c(0, 40), main="Low water")
z = groupdata$treatment=="high_nut"
plot(gen2_length_longest_stem~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5, ylim = c(0, 40), main="High nutrients")
z = groupdata$treatment=="low_nut"
plot(gen2_length_longest_stem~BIO12, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5, ylim = c(0, 40), main="Low nutrients")
legend("topright", c("10", "50", "150", "200", "250"), col=c(10, 50, 150,200,250), pch=15, title="Fitness (seeds)")

# Stem length versus ffdate, by treatment with fitness color scheme
par(mfrow=c(2,2))
palette(rain.colors(250)); 
z = groupdata$treatment=="high_water"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="High water")
z = groupdata$treatment=="low_water"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="Low water")
z = groupdata$treatment=="high_nut"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="High nutrients")
z = groupdata$treatment=="low_nut"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="Low nutrients")
legend("topright", c("10", "50", "150", "250"), col=c(10, 50, 150,250), pch=15, title="Fitness (seeds)", cex=1.25)


# SAME (Stem length versus ffdate) but with Precip color scheme
par(mfrow=c(2,2))
palette(rain.colors(1000)); 
z = groupdata$treatment=="high_water"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="High water")
z = groupdata$treatment=="low_water"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="Low water")
z = groupdata$treatment=="high_nut"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="High nutrients")
z = groupdata$treatment=="low_nut"
plot(gen2_length_longest_stem~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 40), main="Low nutrients")
legend("topright", c("10", "250", "500", "1000"), col=c(10, 250,500,1000), pch=15, title="Precip (mm)", cex=1.25)




# Fitness versus ffdate, but with Precip color scheme
par(mfrow=c(2,2))
palette(rain.colors(1000)); 
z = groupdata$treatment=="high_water"
plot(gen2_no_seeds~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), main="High water")
z = groupdata$treatment=="low_water"
plot(gen2_no_seeds~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), main="Low water")
z = groupdata$treatment=="high_nut"
plot(gen2_no_seeds~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), main="High nutrients")
z = groupdata$treatment=="low_nut"
plot(gen2_no_seeds~days_to_flower, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), main="Low nutrients")
legend("topright", c("10", "250", "500", "1000"), col=c(10, 250,500,1000), pch=15, title="Precip (mm)", cex=1.25)

# All treatments together
plot(gen2_no_seeds~days_to_flower, data=groupdata,col = trmtcols, pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), ylab="population mean fitness", xlab="population mean days to flower", cex.axis=1.2, cex.lab=1.5)
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="high_nut",])), lwd=2, col="green")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="high_water",])), lwd=2, col="blue")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="low_nut",])), lwd=2, col="orange")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="low_water",])), lwd=2, col="yellow")
legend("topright", c("High nutrients", "High water", "Low nutrients", "Low water"), col=c( "green", "blue", "orange", "yellow"), pch=rep(16,4), cex=1.25)


# All treatments together
palette("default")
plot(gen2_no_seeds~days_to_flower, data=groupdata,col = (groupdata$BIO12<360)+1, pch=16, cex=2, xlim=c(50,150), ylim = c(0, 100), ylab="population mean fitness", xlab="population mean days to flower", cex.axis=1.2, cex.lab=1.5)
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="high_nut",])), lwd=2, col="green")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="high_water",])), lwd=2, col="blue")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="low_nut",])), lwd=2, col="orange")
abline(coef(lm(gen2_no_seeds~days_to_flower, data=groupdata[groupdata$treatment=="low_water",])), lwd=2, col="yellow")
legend("topright", c("10", "250", "500", "1000"), col=c(10, 250,500,1000), pch=15, title="Precip (mm)", cex=1.25)


# All data
plot(log(gen2_no_seeds+1)~days_to_flower, data=d.g,col =trmtcols.all, pch=16, cex=1.1,  ylab="Seeds produced", xlab="Days to flower", cex.axis=1.2, cex.lab=1.5)
legend("topright", c("High nutrients", "High water", "Low nutrients", "Low water"), col=c( "green", "blue", "orange", "black"), pch=rep(16,4), cex=1.25)

### HERE CHECK FOR QUADRATIC FITNESS RESPONSE
d.g$days_to_flower2 = d.g$days_to_flower^2
d.g$days_to_flower3 = d.g$days_to_flower^3
days = seq(40, 160, by=1)

summary(lm(log(gen2_no_seeds+1)~days_to_flower + days_to_flower2, data=d.g[d.g$treatment=="low_nut",]))
plot(days, 2.95+ 0.03*days - 0.0003*days^2)

summary(lm(log(gen2_no_seeds+1)~days_to_flower + days_to_flower2+days_to_flower3, data=d.g[d.g$treatment=="high_nut",]))
plot(days, 6.32 - 0.003*days - 0.0003*days^2)

summary(lm(log(gen2_no_seeds+1)~days_to_flower + days_to_flower2, data=d.g[d.g$treatment=="high_water",]))
plot(days, 3.02+0.056*days - 0.00052*days^2)

summary(lm(log(gen2_no_seeds+1)~days_to_flower + days_to_flower2+days_to_flower3, data=d.g[d.g$treatment=="low_water",]))
plot(days,2 + 0.037*days - 0.00032*days^2)

# There appears to be an optimum flowring data in all but the high-nutrient treatment
lines(days, 2.95+0.03*days - 0.0003*days^2, ylim=c(0, 6),type="l", lwd=2, col="black")
lines(days, 3.02+0.056*days - 0.00052*days^2, lwd=2, col="blue")
lines(days, 2+0.037*days - 0.00032*days^2, lwd=2, col="orange")
lines(days, 6.32-0.003*days - 0.0003*days^2, lwd=2, col="green")


m = gam(log(gen2_no_seeds+1)~s(days_to_flower, 2), data=d.g[d.g$treatment=="high_water",], na.action=na.omit)
plot(m, se=T)
summary(m)
m
# So there is some association with Precip and fitness (and long stem length) but not super strong. Other environmental variables? 
summary(glm(gen2_no_seeds~days_to_flower+BIO12, groupdata, family="poisson"))
summary(lm(log(gen2_no_seeds)~days_to_flower+BIO12*Erodium_total + Erodium_total+treatment, groupdata))



plot(days_to_flower~gen2_length_longest_stem, data=groupdata[z,],col = groupdata$gen2_no_seeds[z], pch=16, cex=1.5)
palette(rain.colors(1000)); plot(days_to_flower~gen2_no_seeds, data=groupdata[z,],col = groupdata$BIO12[z], pch=16, cex=1.5)

plot(gen2_length_longest_stem~precipCV, groupdata, col=trmtcols, pch=16)
summary(lm(gen2_length_longest_stem~BIO12+precipCV, groupdata))

plot(days_to_flower~precipCV, groupdata, col=trmtcols, pch=16)
summary(lm(days_to_flower~BIO12+precipCV, groupdata))


plot(waterdiff_internodes~precipCV, groupplastdata, pch=16)
summary(lm(gen2_length_longest_stem~BIO12+precipCV, groupdata))


# NExt to check will be plasticity in this growth pattern response, and whether that's related to source site characteristics. 

# Make a data set with one plasticity value for each pair of treatments (high vs low water, high vs low nutrients)

groupplastdata = groupdata[groupdata$treatment %in% c("high_water"),]

numeric.cols = c(1:3, 13:98) # which columns are numeric and thus should be differenced
groups = unique(d.g$Site)
ngroups = length(groups)
treatments = unique(d.g$treatment)[1:4]
waterdiff_days_to_flower = rep(NA, ngroups)
nutdiff_days_to_flower = rep(NA, ngroups)
for (i in 1:ngroups) {
	groupindex = groupdata$Site == groups[i]
	waterdiff_internodes[i] = groupdata$gen2_no_internodes[groupindex & groupdata$treatment=="high_water"] - groupdata$gen2_no_internodes[groupindex & groupdata$treatment=="low_water"]
	nutdiff_internodes[i] = groupdata$gen2_no_internodes[groupindex & groupdata$treatment=="high_nut"] - groupdata$gen2_no_internodes[groupindex & groupdata$treatment=="low_nut"]
		waterdiff_days_to_flower[i] = groupdata$days_to_flower[groupindex & groupdata$treatment=="high_water"] - groupdata$days_to_flower[groupindex & groupdata$treatment=="low_water"]
	nutdiff_days_to_flower[i] = groupdata$days_to_flower[groupindex & groupdata$treatment=="high_nut"] - groupdata$days_to_flower[groupindex & groupdata$treatment=="low_nut"]
		waterdiff_elong[i] = groupdata$length_per_internode[groupindex & groupdata$treatment=="high_water"] - groupdata$length_per_internode[groupindex & groupdata$treatment=="low_water"]
	nutdiff_elong[i] = groupdata$length_per_internode[groupindex & groupdata$treatment=="high_nut"] - groupdata$length_per_internode[groupindex & groupdata$treatment=="low_nut"]
}
	groupplastdata = cbind(groupplastdata, waterdiff_internodes, waterdiff_elong, waterdiff_days_to_flower, nutdiff_internodes, nutdiff_elong, nutdiff_days_to_flower)


# Plot these plasticity values against environmental factors

plot(waterdiff_internodes~BIO12, groupplastdata, col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)")
points(nutdiff_internodes~BIO12, groupplastdata, col = "orange", pch=16, cex=1.5)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))
abline(c(0,0), lwd=2, lty=2, col=gray(0.5))


par(mfrow=c(1,2))
plot(waterdiff_internodes~BIO12, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="California", ylab="Plasticity in number of internodes", ylim=c(0,5))
points(nutdiff_internodes~BIO12, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=1.5)
plot(waterdiff_internodes~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="Chile", ylab="Plasticity in number of internodes", ylim=c(0,5))
points(nutdiff_internodes~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=1.5)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))



plot(waterdiff_days_to_flower~BIO12, groupplastdata, col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", ylim = c(-40, 10))
points(nutdiff_days_to_flower~BIO12, groupplastdata, col = "orange", pch=16, cex=1.5)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))
abline(c(0,0), lwd=2, lty=2, col=gray(0.5))



par(mfrow=c(1,2))
plot(waterdiff_days_to_flower~BIO12, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="Mean annual precip (mm)", main="California", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~BIO12, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=2)
plot(waterdiff_days_to_flower~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="Mean annual precip (mm)", main="Chile", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))



# With respect to source-site total cover -- some trend toward more plasticity in first flower date associated with higher total cover. But not significant in this analysis. 
par(mfrow=c(1,2))
plot(waterdiff_days_to_flower~total_cover, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="Erodium cover", main="California", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~total_cover, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=2)
plot(waterdiff_days_to_flower~total_cover, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="Erodium cover", main="Chile", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~total_cover, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))

summary(lm(nutdiff_days_to_flower~total_cover, groupplastdata, subset=groupplastdata$region=="California"))

# With respect to source-site mean Phosphorus -- hint of more plasticity with respect to nutrients from low-P sites.
par(mfrow=c(1,2))
plot(waterdiff_days_to_flower~P1, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="Phosphorus (ppm)", main="California", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~P1, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=2)
plot(waterdiff_days_to_flower~P1, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="Phosphorus (ppm)", main="Chile", ylab="Plasticity in flowering date (days)")
points(nutdiff_days_to_flower~P1, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))

# Check plasticity vs within-site variation in P1
CVP1 = rep(NA, nsites)
for (i in 1:nsites) {
	z = d.g$Site==sites[i]
	CVP1[i] = var(d.g$P1[z], na.rm=T)/mean(d.g$P1[z], na.rm=T)
}
z = match(groupplastdata$Site, sites)
groupplastdata$CVP1 = CVP1[z]

# With respect to within-source-site CV Phosphorus -- no pattern with flowering date. Nothing clear with internodes either. 
par(mfrow=c(1,2))
plot(waterdiff_internodes~CVP1, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="Phosphorus (ppm)", main="California", ylab="Plasticity in internodes")
points(nutdiff_internodes~CVP1, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=2)
plot(waterdiff_internodes~CVP1, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="Phosphorus (ppm)", main="Chile", ylab="Plasticity in internodes")
points(nutdiff_internodes~CVP1, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))



# Check plasticity vs within-site variation in total cover
CVcover = rep(NA, nsites)
for (i in 1:nsites) {
	z = d.g$Site==sites[i]
	CVcover[i] = var(d.g$total_cover[z], na.rm=T)/mean(d.g$P1[z], na.rm=T)
}
z = match(groupplastdata$Site, sites)
groupplastdata$CVcover = CVcover[z]

# With respect to within-source-site CV total cover -- some trend toward lower plasticity in internodes with higher CV cover. No clear trend for flower date. 
par(mfrow=c(1,2))
plot(waterdiff_days_to_flower~CVcover, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=2, xlab="CV cover", main="California", ylab="Plasticity in internodes")
points(nutdiff_days_to_flower~CVcover, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=2)
plot(waterdiff_days_to_flower~CVcover, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=2, xlab="CV cover", main="Chile", ylab="Plasticity in internodes")
points(nutdiff_days_to_flower~CVcover, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=2)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))









plot(waterdiff_elong~BIO12, groupplastdata, col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", ylab="Plasticity in internode length")
points(nutdiff_elong~BIO12, groupplastdata, col = "orange", pch=16, cex=1.5)
legend("topright", c("plasticity to water", "plasticity to nutrients"), col=c("slateblue", "orange"), pch=rep(16,2))
abline(c(0,0), lwd=2, lty=2, col=gray(0.5))

par(mfrow=c(1,2))
plot(waterdiff_elong~BIO12, groupplastdata[groupplastdata$region=="California",], col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="California")
points(nutdiff_elong~BIO12, groupplastdata[groupplastdata$region=="California",], col = "orange", pch=16, cex=1.5)
plot(waterdiff_elong~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "slateblue",pch=16, cex=1.5, xlab="Mean annual precip (mm)", main="Chile")
points(nutdiff_elong~BIO12, groupplastdata[groupplastdata$region=="Chile",], col = "orange", pch=16, cex=1.5)
#legend("bottomright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))

# Check corrected days to flower data

summary(lm(days_to_flower~gen2_length_longest_stem+BIO12+P1+total_cover, d.g))

summary(lm(days_to_flower~gen2_length_longest_stem+BIO12+P1+total_cover, groupdata))
cor(groupdata$days_to_flower, groupdata$gen2_length_longest_stem)
plot(days_to_flower~gen2_length_longest_stem, groupdata, pch=16, col=trmtcols)
plot(days_to_flower~gen2_length_longest_stem, groupdata, pch=16, col=groupdata$BIO12)


plot(days_to_flower~BIO12, groupdata, col = trmtcols,pch=16, cex=1.5, xlab="Mean annual precip (mm)")
legend("topright", c("High water", "Low water", "High nutrients", "Low nutrients"), col=c("blue", "yellow", "green", "orange"), pch=rep(16,4))
