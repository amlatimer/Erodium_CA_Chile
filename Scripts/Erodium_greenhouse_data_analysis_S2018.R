# Erodium cicutarium greenhouse experiment analysis for re-submission of paper 
# Andrew Latimer
# May 2018

# Note this uses absolute plasticity values

###########################
# Load data and libraries
###########################

library(lattice); library(fields); library(MASS); library(reldist); library(lme4); library(lmerTest); library(MuMIn); library(nlme); library(lmtest); library(ggplot2); library(lars); library(glmnet); library(hier.part); library(corrplot); library(plyr); library(dplyr)

source("Ecic_analysis_functions.R")

d = read.csv("../Data/All field data combined.csv")
d.g = read.csv("../Data/second generation dataset.csv")
popdata = read.csv("../Data/allpopsinfo.csv")
sitedata = read.csv("../Data/Sitedata.csv")
#fielddata = read.csv("../Data/All field data combined.csv")

# read in and add pptcv information from Brooke and Ernesto -- measure of temporal heterogeneity
pptcv = read.csv("../Data/pptcv_by_site.csv")
d.g = merge(d.g, pptcv, by="Site")
head(d.g)


##############################################
# Check key greenhouse data variables for errors and outliers
#############################################
pairs(d.g[,c("days_to_flower", "gen2_length_longest_stem", "SLA", "Total.Area..mm2.", "gen2_total_seeds")])
stem(d.g$gen2_length_longest_stem)
# There appears to be one outlier in stem length (1 value > 100 while other values all below 70)
# remove this one huge outlier, appears to represent a typo
d.g = d.g[-which.max(d.g$gen2_length_longest_stem),]

# check SLA values 
# convert to mm2/mg units
hist(d.g$SLA/1000) 
d.g$SLA = d.g$SLA / 1000

# there should be no zero stem lengths, so treat these values as missing data and replace with NA
d.g$gen2_length_longest_stem[d.g$gen2_length_longest_stem==0] = NA


###########################################################
# CHECKING BODEGA BAY
################################################
# Bodega Bay cove is an outlier for mortality overall and in the low-water treatment especially: only 33/80 plants survived to set seed. This leads to extreme imbalance across treatments and populations. Therefore we removed this population from the data set for all analyses.
d.g.complete = filter(d.g, gen2_no_seeds>0)
table(d.g.complete$Site, d.g.complete$treatment)
d.g = d.g[d.g$Site != "Bodega Bay cove", ]
d.g = droplevels(d.g)

###############################################################
# Analysis by linear mixed models of greenhouse responses: Trait means
# Addresses Q1 in manuscript: 
# Are there differences in traits means among regions, populations and treatments? 
###############################################################

###### First, compute family means for each treatment 

d.g.fam_numbers = aggregate(d.g[,13:ncol(d.g)], by=list(d.g$treatment, d.g$field_family_ID), FUN=mean, na.rm=T)
d.g.fam_labels = aggregate(d.g[,1:12], by=list(d.g$treatment, d.g$field_family_ID), FUN=f<-function(x){return(x[1])})
d.g.fam = cbind(d.g.fam_labels, d.g.fam_numbers)
head(d.g.fam)

# There are a few families with leaf area == 0 - since this doesn't make sense, switch these to NAs
d.g.fam$Total.Area..mm2.[which(d.g.fam$Total.Area..mm2.==0)] = NA
table(d.g.fam$Site)

#######################################
#### Create soil principal components
#######################################
soildata = d.g.fam[,c("P1", "OM", "PH", "K", "MG", "CA", "CEC")]
soildata=scale(soildata)
soilpc_rows = complete.cases(soildata)
soil_pc = prcomp(soildata[soilpc_rows,])
biplot(soil_pc)
soil_pc
summary(soil_pc)
# First PC reflects mainly soil organic matter and CEC as well as MG and CA content. No single PC reflects majority of variation. But because first PC is broad measure of nutrient profile, seems OK to use it as the local soil variability measure.
# Add this first PC into the data set 
d.g.fam$soil_pc1 = NA
d.g.fam$soil_pc2 = NA
d.g.fam$soil_pc1[soilpc_rows] = soil_pc$x[,1]
d.g.fam$soil_pc2[soilpc_rows] = soil_pc$x[,2]

############################
# Calculate family-specific plasticity values

families = unique(d.g$field_family_ID)
n_families=length(families)
famplast.water = d.g[match(families, d.g$field_family_ID),]
famplast.nut = d.g[match(families, d.g$field_family_ID),]

for (i in 1:n_families) {
  
  famplast.water$d2f_mean[i] = mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment =="high_water"], na.rm=T)
  famplast.nut$d2f_mean[i] = mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment =="high_nut"], na.rm=T)
  
  famplast.water$sl_mean[i] = mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment =="high_water"], na.rm=T)
  famplast.nut$sl_mean[i] = mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment =="high_nut"], na.rm=T)
  
  famplast.water$area_mean[i] = mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment =="high_water"], na.rm=T)
  famplast.nut$area_mean[i] = mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment =="high_nut"], na.rm=T)
  
  famplast.water$SLA_mean[i] = mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment =="high_water"], na.rm=T)
  famplast.nut$SLA_mean[i] = mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment =="high_nut"], na.rm=T)
  
  famplast.water$d2f_plast[i] = abs((mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment=="low_water"], na.rm=T) - mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment=="high_water"], na.rm=T)) / famplast.water$d2f_mean[i])
  famplast.nut$d2f_plast[i] = abs((mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment=="low_nut"], na.rm=T) - mean(d.g$days_to_flower[d.g$field_family_ID==families[i] & d.g$treatment=="high_nut"], na.rm=T)) / famplast.nut$d2f_mean[i])
  
  famplast.water$sl_plast[i] = abs((mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment=="low_water"], na.rm=T) - mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment=="high_water"], na.rm=T)) / famplast.water$sl_mean[i])
  famplast.nut$sl_plast[i] = abs((mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment=="low_nut"], na.rm=T) - mean(d.g$gen2_length_longest_stem[d.g$field_family_ID==families[i] & d.g$treatment=="high_nut"], na.rm=T)) / famplast.nut$sl_mean[i])
  
  famplast.water$area_plast[i] = abs((mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment=="low_water"], na.rm=T) - mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment=="high_water"], na.rm=T)) / famplast.water$area_mean[i])
  famplast.nut$area_plast[i] = abs((mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment=="low_nut"], na.rm=T) - mean(d.g$Total.Area..mm2.[d.g$field_family_ID==families[i] & d.g$treatment=="high_nut"], na.rm=T)) / famplast.nut$area_mean[i])
  
  famplast.water$SLA_plast[i] = abs((mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment=="low_water"], na.rm=T) - mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment=="high_water"], na.rm=T)) / famplast.water$SLA_mean[i])
  famplast.nut$SLA_plast[i] = abs((mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment=="low_nut"], na.rm=T) - mean(d.g$SLA[d.g$field_family_ID==families[i] & d.g$treatment=="high_nut"], na.rm=T)) / famplast.nut$SLA_mean[i])
  
}

#########################################################################
# To streamline analysis, combine all data sets into a single data frame
#   with one row for each family,  
#   one column for each trait mean value across all treatments, 
#   and one column for each trait plasticity value. 
#########################################################################

# Assemble trait means for each family across all treatments
# select only the columns needed for analysis
traits = c("days_to_flower", "gen2_length_longest_stem", "SLA", "Total.Area..mm2.")
cols_to_use = c("Site", "region", "field_family_ID", "treatment", traits)
dsub = d.g.fam[,cols_to_use]
d.traits = ddply(dsub, .(field_family_ID), summarize, Site=unique(Site), region=unique(region), days_to_flower_mean=mean(days_to_flower), stem_length_mean=mean(gen2_length_longest_stem), SLA_mean=mean(SLA), leaf_area_mean=mean(Total.Area..mm2.))
head(d.traits)

# Assemble plasticity values 
traits2 = c("d2f_plast", "sl_plast", "SLA_plast", "area_plast")
cols_to_use = c("field_family_ID", traits2)
plastsub_water = famplast.water[,cols_to_use]
names(plastsub_water)[2:5] = c("d2f_plast_water", "sl_plast_water", "SLA_plast_water", "area_plast_water")
plastsub_nut = famplast.nut[,cols_to_use]
names(plastsub_nut)[2:5] = c("d2f_plast_nut", "sl_plast_nut", "SLA_plast_nut", "area_plast_nut")

# Put all the data together
d.traits.all = merge(d.traits, plastsub_water, by="field_family_ID")
d.traits.all = merge(d.traits.all, plastsub_nut, by="field_family_ID")
head(d.traits.all)

# Use mixed models to test for differences between treatments, populations, and regions for trait means and plasticity values. Treatments are fixed effect, Site within region are nested random effects. 
# Model testing done using lmerTest() 

# Test example -- Fit mixed models
m1.d2f = lmer(days_to_flower~treatment+(1|region), data=d.g.fam[!is.na(d.g.fam$days_to_flower),], REML=FALSE)
m2.d2f = lmer(days_to_flower~treatment+(1|region/Site), data=d.g.fam[!is.na(d.g.fam$days_to_flower),], REML=FALSE)
# Test significance of region random effects
rand(m1.d2f)
# Test significance of Site random effects
anova(m1.d2f, m2.d2f)
# Test significance of treatment fixed effect
anova(m2.d2f) 


##############################################
# Apply model testing framework to all trait means and to plasticity values
##############################################

test_trait_plast <- function(d, traits) { 
  # Fit mixed-effects models using lmer() then do model comparison using lmerTest. 
  # Full model is of form: trait ~ treatment + (1|region/Site)
  # Tests significance of region random effect, site random effect nested within region, and treatment fixed effect. 
  # Arguments d is a data frame containing columns including family mean trait values for each treatment, a treatment indicator column, a Site (or population) indicator column, and a region indicator column. 
  # Argument traits is a vector of column names of the traits to analyze
  n.traits = length(traits)
  trait_labels = chisq_region = chisq_site = pval_region = pval_site = {}

  for (i in 1:n.traits) {
    y = scale(d[,traits[i]])
    dtemp = data.frame(trait=y, region=d$region, Site=d$Site)
    # remove any rows containing NAs in any of the variables
    dtemp = dtemp[complete.cases(dtemp),]
    # fit models
    fit_full = lmer(trait ~ (1 | region/Site), data=dtemp, REML=FALSE)
    fit_reduced = lmer(trait~ (1|region), data=dtemp, REML=FALSE)
    # model tests
    lrtest_region = rand(fit_reduced)
    lrtest_site = anova(fit_reduced, fit_full)
    # pull out and store results 
    trait_labels = c(trait_labels, traits[i])
    chisq_region = c(chisq_region, lrtest_region$rand.table$Chi.sq)
    pval_region = c(pval_region, lrtest_region$rand.table$p.value)
    chisq_site = c(chisq_site, lrtest_site$Chisq[2])
    pval_site = c(pval_site, lrtest_site$`Pr(>Chisq)`[2])
  }
  
  out_table = data.frame(trait=trait_labels, chisq_region, pval_region, chisq_site, pval_site)
  # adjusted p-values
  n.models = length(pval_region)
  pval_vec = c(pval_region, pval_site)
  qval_vec = p.adjust(pval_vec, method="fdr")
  out_table$pval_region_adj = qval_vec[1:n.models]
  out_table$pval_site_adj = qval_vec[(n.models+1):(2*n.models)]
  # reorder output so it can be pasted directly into table for the paper
  out_table = out_table[,c("trait", "chisq_region", "pval_region_adj", "chisq_site", "pval_site_adj")]
  
  return(out_table)
}

test_trait_means <- function(d, traits) { 
  # Fit mixed-effects models using lmer() then do model comparison using lmerTest. 
  # Full model is of form: trait ~ (1|region/Site)
  # Tests significance of region random effect and site random effect nested within region.
  # Arguments d is a data frame containing columns including family mean plasticity values for each trait, a Site (or population) indicator column, and a region indicator column. 
  # Argument traits is a vector of column names of the plasticity traits to analyze
  n.traits = length(traits)
  trait_labels = chisq_region = chisq_site = F_treatment =  pval_region = pval_site = pval_treatment = {}
  
  for (i in 1:n.traits) {
    y = scale(d[,traits[i]])
    dtemp = data.frame(trait=y, treatment=d$treatment, region=d$region, Site=d$Site)
    # remove any rows containing NAs in any of the variables
    dtemp = dtemp[complete.cases(dtemp),]
    # fit models
    fit_full = lmer(trait ~ treatment + (1 | region/Site), data=dtemp, REML=FALSE)
    fit_reduced = lmer(trait~treatment + (1|region), data=dtemp, REML=FALSE)
    # model tests
    anova_treatment = anova(fit_full)
    lrtest_region = rand(fit_reduced)
    lrtest_site = anova(fit_reduced, fit_full)
    # pull out and store results 
    trait_labels = c(trait_labels, traits[i])
    chisq_region = c(chisq_region, lrtest_region$rand.table$Chi.sq)
    pval_region = c(pval_region, lrtest_region$rand.table$p.value)
    chisq_site = c(chisq_site, lrtest_site$Chisq[2])
    pval_site = c(pval_site, lrtest_site$`Pr(>Chisq)`[2])
    F_treatment = c(F_treatment, anova_treatment$F.value)
    pval_treatment = c(pval_treatment, anova_treatment$`Pr(>F)`)
  }
  
  out_table = data.frame(trait=trait_labels, chisq_region, pval_region, chisq_site, pval_site, F_treatment, pval_treatment)
  # adjusted p-values
  n.models = length(pval_treatment)
  pval_vec = c(pval_region, pval_site, pval_treatment)
  qval_vec = p.adjust(pval_vec, method="fdr")
  out_table$pval_region_adj = qval_vec[1:n.models]
  out_table$pval_site_adj = qval_vec[(n.models+1):(2*n.models)]
  out_table$pval_treatment_adj = qval_vec[(2*n.models+1):(3*n.models)]
  # reorder output so it can be pasted directly into table for the paper
  out_table = out_table[,c("trait", "chisq_region", "pval_region_adj", "chisq_site", "pval_site_adj", "F_treatment", "pval_treatment_adj")]
  return(out_table)
}

### Do the analysis for trait means
traits_means = c("days_to_flower", "gen2_length_longest_stem", "SLA", "Total.Area..mm2.")
trait_results_means = test_trait_means(d.g.fam, traits_means)
trait_results_means

    ### Do the analysis for trait plasticity values
traits_plast = names(d.traits.all)[grep("plast",names(d.traits.all))]
trait_results_plast = test_trait_plast(d.traits.all, traits_plast)
trait_results_plast

# Put together into table format for paper
trait_results_plast$F_treatment = trait_results_plast$pval_treatment_adj = NA
Q1_results = rbind(trait_results_means, trait_results_plast)

write.csv(Q1_results, "../Results/Q1_results.csv")

########################################################
# Plot Figure 2
########################################################

### Figure showing variance explained by treatment, site, interaction
par(mfrow=c(1,1))
m.d2f <- lm(days_to_flower ~ treatment*Site, data=d.g.fam)
m.sl <- lm(gen2_length_longest_stem ~ treatment*Site, data=d.g.fam)
m.sla <- lm(SLA ~ treatment*Site, data=d.g.fam)
m.area <- lm(Total.Area..mm2. ~ treatment*Site, data=d.g.fam)

a = list(m.d2f, m.sl, m.sla, m.area)
a.names = c("Days to \nflower", "Stem\nlength", "SLA", "Leaf\narea")
lapply(a, FUN=anova)
aov.stackbar(a, a.names) # Figure showing percent variance explained for each response variable



##############################
# Make interaction plot (Figure 3)

pdf("../Plots/Figure_3_reaction_norms_water.pdf")
# Optionally color the lines by source-site precipitation
#grayscale <- colorRampPalette(c("light gray",  "black"))
#color_var = sitedata$had_pptmean
#max_col = 100
#cols = (color_var - min(color_var)) * (99/ (max(color_var) - min(color_var))) + 1
#palette(grayscale(max_col))

# Alternatively color the lines by source region
palette("default")

cols = sapply(as.character(sitedata$region), switch, "Chile" = gray(0.7), "California" = gray(0.3))
axis.scale = 1.6
par(mfrow=c(2,2), mar=c(5,5,5,5))
# days to flower, plasticity to water 
plotdata = d.g.fam[d.g.fam$treatment %in% c("high_water","low_water"),]
plotdata$treatment = droplevels(plotdata$treatment)
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$days_to_flower, type="b", xlab="Treatment", ylab="Days to flowering", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# stem length, plasticity to water
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$gen2_length_longest_stem, type="b", xlab="Treatment", ylab="Length of longest stem (cm)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# SLA, plasticity to water
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$SLA, type="b", xlab="Treatment", ylab="SLA (mm2/mg)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# leaf area, plasticity to water
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$Total.Area..mm2., type="b", xlab="Treatment", ylab="Leaf area (mm2)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)

dev.off()

pdf("../Plots/Figure_3_reaction_norms_nutrients.pdf")

# days to flower, plasticity to nutrients 
if (names(d.g.fam)[1] == "Group.1") d.g.fam = d.g.fam[,-c(1,2)]
plotdata = filter(d.g.fam, treatment %in% c("high_nut","low_nut"))
plotdata$treatment = droplevels(plotdata$treatment)
par(mfrow=c(2,2))
# days to flower, plasticity to nutrients 
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$days_to_flower, type="b", xlab="Treatment", ylab="Days to flower", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# stem length, plasticity to nutrients
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$gen2_length_longest_stem, type="b", xlab="Treatment", ylab="Length of longest stem (cm)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# SLA, plasticity to nutrients
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$SLA, type="b", xlab="Treatment", ylab="SLA (mm2/mg)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)
# leaf area, plasticity to nutrients
interaction.plot(x.factor=plotdata$treatment, trace.factor=plotdata$Site, response=plotdata$Total.Area..mm2., type="b", xlab="Treatment", ylab="Leaf area (mm2)", cex.axis=axis.scale, cex.lab=axis.scale, legend=FALSE, function(x) mean(x, na.rm = TRUE), pch=16, lty=1, lwd=2, col=cols)

dev.off()



#################################################
# Question 2 -- assocations between population-level trait means and plasticity levels and source-site environments
#################################################


#######################################################################
# load data sets that include the hadcrut pptcv data 
haddata = read.csv("../data/hadley_ppt_data_by_site.csv")
worldclimdata = read.csv("../data/popworldclimdata.csv")


# Summarize trait means and plasticity values to site-level means

traits.all.site = ddply(d.traits.all, .(Site), summarize, region=unique(region), days_to_flower_mean=mean(days_to_flower_mean, na.rm=T), stem_length_mean=mean(stem_length_mean, na.rm=T), SLA_mean=mean(SLA_mean, na.rm=T), leaf_area_mean=mean(leaf_area_mean, na.rm=T), d2f_plast_water=mean(d2f_plast_water, na.rm=T), sl_plast_water=mean(sl_plast_water, na.rm=T), SLA_plast_water=mean(SLA_plast_water, na.rm=T), area_plast_water=mean(area_plast_water, na.rm=T), d2f_plast_nut=mean(d2f_plast_nut, na.rm=T), sl_plast_nut=mean(sl_plast_nut, na.rm=T), SLA_plast_nut=mean(SLA_plast_nut, na.rm=T), area_plast_nut=mean(area_plast_nut, na.rm=T))
head(traits.all.site)

# create data sets for trait means that are specific to a particular treatment
if(names(d.g.fam)[1]=="Group.1") d.g.fam = d.g.fam[,-c(1,2)] # remove extra grouping columns left from aggregate()
sitedata_lowater = filter(d.g.fam, treatment=="low_water") %>%  ddply(.(Site), summarize, region=unique(region), days_to_flower_mean=mean(days_to_flower, na.rm=T), stem_length_mean=mean(gen2_length_longest_stem, na.rm=T), SLA_mean=mean(SLA, na.rm=T), leaf_area_mean=mean(Total.Area..mm2., na.rm=T))
head(sitedata_lowater)

sitedata_highwater = filter(d.g.fam, treatment=="high_water") %>%  ddply(.(Site), summarize, region=unique(region), days_to_flower_mean=mean(days_to_flower, na.rm=T), stem_length_mean=mean(gen2_length_longest_stem, na.rm=T), SLA_mean=mean(SLA, na.rm=T), leaf_area_mean=mean(Total.Area..mm2., na.rm=T))

# Summarize the plant-level environmental variables to their mean and CV by site. 
#   This includes the soil and plant community cover data. 
cv <- function(x) {return(sd(x, na.rm=T)/mean(x, na.rm=T))}
field_data_site = ddply(d.g.fam, .(Site), summarize, soil_pc1_mean=mean(soil_pc1, na.rm=T), soil_pc1_cv=cv(soil_pc1), soil_pc2_mean=mean(soil_pc2, na.rm=T), soil_pc2_cv=cv(soil_pc2),total_cover_CV=cv(total_cover), total_cover=mean(total_cover, na.rm=T))
head(field_data_site)

# merge with Hadley data
traits.all.site = merge(traits.all.site, haddata, by="Site")
head(traits.all.site)
sitedata_lowater = merge(sitedata_lowater, haddata, by="Site")
sitedata_highwater = merge(sitedata_highwater, haddata, by="Site")
head(sitedata_highwater)

# merge with Worldclim data
names(worldclimdata)[1] = "Site"
traits.all.site = merge(traits.all.site, worldclimdata[,c("Site", "BIO1", "BIO12", "BIO15")], by="Site")
head(traits.all.site)
sitedata_lowater = merge(sitedata_lowater, worldclimdata[,c("Site", "BIO1", "BIO12", "BIO15")], by="Site") 
sitedata_highwater = merge(sitedata_highwater, worldclimdata[,c("Site", "BIO1", "BIO12", "BIO15")], by="Site")
head(sitedata_highwater)

# merge with family-level field data (soils, cover)
traits.all.site = merge(traits.all.site, field_data_site, by="Site")
head(traits.all.site)
sitedata_lowater = merge(sitedata_lowater, field_data_site, by="Site")
sitedata_highwater = merge(sitedata_highwater, field_data_site, by="Site")
head(sitedata_highwater)


###################################################################
# Test one explanatory variable at a time, and check for interaction with region 
###################################################################
allx = c("BIO1", "BIO12", "BIO15", "had_pptcv", "soil_pc1_mean", "soil_pc1_cv", "total_cover", "total_cover_CV")


### Q HERE DO WE NEED TO TEST FOR INTERACTION WITH REGION ONE VARIABLE AT A TIME, BEFORE TESTING THE COMBINED MODELS? 



#########################################################################
# global model search and ranking via MuMIn

# For results section, make plots of importance values, and tables of coefficients for the "best" model. 
# All explanatory variables 
# Note I tried analyzing while also including the second soil principal component and it was marginally important and didn't change the main results. So for simplicity, I am using only the first soil PC. 

# Explanatory variables to include in the analysis
allx = c("BIO1", "BIO12", "BIO15", "had_pptcv", "soil_pc1_mean", "soil_pc1_cv", "total_cover", "total_cover_CV")
# Names of the traits for labeling plots
xnames = c("MAT", "MAP", "PPTSEAS", "PPTCV", "soil_nut", "CV_soil_nut", "total_cover", "CV_total_cover")

# Define the site-level data set to work on 
sitedata=traits.all.site

# First transform some of the response variables
sitedata$log_stemlength = log(sitedata$stem_length_mean)
sitedata$log_sla = log(sitedata$SLA_mean)
sitedata$log_area = log(sitedata$leaf_area_mean)

### Set na.action for MuMIn
options(na.action = "na.fail")

# Set up the cluster for running pdredge()
library(snow)
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 4), type = clusterType))

# Make a plot of importance values 

pdf("../Plots/Importance_trait_means.pdf")

# Group plots into sets of 4
par(mfrow=c(2,2), mar=c(7,4,2,2)+0.1)
axlab = 1.2
axtitle = 1.1

# days to flowering
testdata = cbind(days_to_flower=sitedata$days_to_flower, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(days_to_flower~., data=testdata)
clusterExport(clust, "testdata")
d2f_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
imp = importance(subset(d2f_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Days to flowering", ylim=c(0,1))
#'Best' modelx
summary(get.models(d2f_ms1, 1)[[1]])
AICc(get.models(d2f_ms1, 2)[[1]]) - AICc(get.models(d2f_ms1, 1)[[1]])
#BIO12

# stem length
testdata = cbind(log_stemlength=sitedata$log_stemlength, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(log_stemlength ~ ., data=testdata)
clusterExport(clust, "testdata")
sl_ms1 = pdredge(fm1, rank = "AICc", cluster=clust)
#plot(sl_ms1)
imp = importance(subset(sl_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Stem length", ylim=c(0,1))
#'Best' model
summary(get.models(sl_ms1, 1)[[1]])
AICc(get.models(sl_ms1, 2)[[1]]) - AICc(get.models(sl_ms1, 1)[[1]])
# BIO12

# SLA
testdata = cbind(log_sla=sitedata$log_sla, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(log_sla ~ ., data=testdata)
clusterExport(clust, "testdata")
SLA_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(SLA_ms1)
imp = importance(subset(SLA_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="SLA", ylim=c(0,1))
#'Best' model
summary(get.models(SLA_ms1, 1)[[1]])
AICc(get.models(SLA_ms1, 2)[[1]]) - AICc(get.models(SLA_ms1, 1)[[1]])
# BIO12, [BIO1]

# Leaf area
testdata = cbind(log_area=sitedata$log_area, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(log_area ~ ., data=testdata)
clusterExport(clust, "testdata")
leafarea_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(leafarea_ms1)
imp = importance(subset(leafarea_ms1,delta <= 5))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Leaf area", ylim=c(0,1))
#'Best' model
summary(get.models(leafarea_ms1, 1)[[1]])
AICc(get.models(leafarea_ms1, 2)[[1]]) - AICc(get.models(leafarea_ms1, 1)[[1]])

dev.off()
# End of plot 


# Make plot of importance values for plasticity to water 
pdf("../Plots/Importance_plasticity_water.pdf")
par(mfrow=c(2,2), mar=c(7,4,2,2)+0.1)
axlab = 1.2
axtitle = 1.1

# plasticity of flowering time to water
testdata = cbind(d2f_plast=sitedata$d2f_plast_water, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(d2f_plast ~., data=testdata)
clusterExport(clust, "testdata")
d2f_wp_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(d2f_wp_ms1) # no variables
imp = importance(subset(d2f_wp_ms1,delta <= 6))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Days to flower", ylim=c(0,1))
#'Best' model
summary(get.models(d2f_wp_ms1, 1)[[1]])
AICc(get.models(d2f_wp_ms1, 2)[[1]]) - AICc(get.models(d2f_wp_ms1, 1)[[1]])


# plasticity of stem length to water
testdata = cbind(sl_plast=sitedata$sl_plast_water, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(sl_plast ~ ., data=testdata)
clusterExport(clust, "testdata")
sl_wp_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(sl_wp_ms1)
imp = importance(subset(sl_wp_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Stem length", ylim=c(0,1))
summary(get.models(sl_wp_ms1, 1)[[1]])
AICc(get.models(sl_wp_ms1, 2)[[1]]) - AICc(get.models(sl_wp_ms1, 1)[[1]])

# plasticity of SLA to water
testdata = cbind(SLA_plast=sitedata$SLA_plast_water, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(SLA_plast ~ ., data=testdata)
clusterExport(clust, "testdata")
SLA_wp_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(SLA_wp_ms1)
imp = importance(subset(SLA_wp_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="SLA", ylim=c(0,1))
summary(get.models(SLA_wp_ms1, 1)[[1]])
AICc(get.models(SLA_wp_ms1, 2)[[1]]) - AICc(get.models(SLA_wp_ms1, 1)[[1]])
# nothing

# plasticity of leaf area to water
testdata = cbind(area_plast=sitedata$area_plast_water, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(area_plast ~ ., data=testdata)
clusterExport(clust, "testdata")
area_wp_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(area_wp_ms1)
imp = importance(subset(area_wp_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Leaf area", ylim=c(0,1))
summary(get.models(area_wp_ms1, 1)[[1]])
AICc(get.models(area_wp_ms1, 2)[[1]]) - AICc(get.models(area_wp_ms1, 1)[[1]])

dev.off()
# End of plot 


# Make plot of importance values for plasticity to nutrients 
pdf("../Plots/Importance_plasticity_nutrients.pdf")
par(mfrow=c(2,2), mar=c(7,4,2,2)+0.1)
axlab = 1.2
axtitle = 1.1

# plasticity of flowering time to nutrients
testdata = cbind(d2f_plast=sitedata$d2f_plast_nut, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(d2f_plast ~ ., data=testdata)
clusterExport(clust, "testdata")
d2f_np_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(d2f_np_ms1)
imp = importance(subset(d2f_np_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=1.3, cex.names = 1.2, main="Days to flower", ylim=c(0,1))
summary(get.models(d2f_np_ms1, 1)[[1]])
AICc(get.models(d2f_np_ms1, 2)[[1]]) - AICc(get.models(d2f_np_ms1, 1)[[1]])

# plasticity of stem length to nutrients
testdata = cbind(sl_plast=sitedata$sl_plast_nut, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] =  xnames
fm1 = lm(log(sl_plast+2) ~ ., data=testdata)
clusterExport(clust, "testdata")
sl_np_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(sl_np_ms1)
imp = importance(subset(sl_np_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Stem length", ylim=c(0,1))
summary(get.models(sl_np_ms1, 1)[[1]])
AICc(get.models(sl_np_ms1, 2)[[1]]) - AICc(get.models(sl_np_ms1, 1)[[1]])

# plasticity of SLA to nutrients
testdata = cbind(SLA_plast=sitedata$SLA_plast_nut, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(SLA_plast ~ ., data=testdata)
#z=step(fm1, direction="backward", k=5)
clusterExport(clust, "testdata")
SLA_np_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(SLA_np_ms1)
imp = importance(subset(SLA_np_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue, names.arg=df$variable, las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="SLA", ylim=c(0,1))
summary(get.models(SLA_np_ms1, 1)[[1]])
AICc(get.models(SLA_np_ms1, 2)[[1]]) - AICc(get.models(SLA_np_ms1, 1)[[1]])


# plasticity of leaf area to nutrients
testdata = cbind(area_plast=sitedata$area_plast_nut, sitedata[,allx])
testdata = testdata[complete.cases(testdata),]
names(testdata)[2:(length(allx)+1)] = xnames
fm1 = lm(area_plast ~., data=testdata)
clusterExport(clust, "testdata")
area_np_ms1 = pdredge(fm1, rank="AICc", cluster=clust)
#plot(area_np_ms1)
imp = importance(subset(area_np_ms1,delta <= 3))
imp = imp[imp>=0.1]
df = data.frame(variable=names(imp), impvalue=matrix(imp))
barplot(df$impvalue[which(df$impvalue>=0.3)], names.arg=df$variable[which(df$impvalue>=0.3)], las=2, ylab="Importance", cex.axis=1.3, cex.lab=axlab, cex.names = axtitle, main="Leaf area", ylim=c(0,1))
summary(get.models(area_np_ms1, 1)[[1]])
AICc(get.models(area_np_ms1, 2)[[1]]) - AICc(get.models(area_np_ms1, 1)[[1]])

dev.off()
# End of plot 

# End of Q2 section
############################################################################


###################################################
# QUESTION 3
# Assess trait correlations and compare them between the two regions 

### Make Figure 4 

d.traits.all = merge(d.traits.all, worldclimdata, by="Site")
names(d.traits.all)[3] = "region"

p1 = ggplot(d.traits.all, aes(x=log(stem_length_mean), y=days_to_flower_mean, color=region)) + guides(color=FALSE) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=12), legend.position="none") +
  geom_point(shape=16) +   
  scale_colour_hue(l=50) + labs(x="log(stem_length cm)", y="days to flowering") +
    geom_smooth(method=lm)    # Add linear regression line 
#  (by default includes 95% confidence region)


p2 = ggplot(d.traits.all, aes(x=log(SLA_mean), y=days_to_flower_mean, color=region)) + guides(color=FALSE) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=12), legend.position="none") +
  geom_point(shape=16) +   
  scale_colour_hue(l=50) + labs(x="log(SLA mm2/mg)", y="days to flowering") +
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)

summary(lm(days_to_flower~log(SLA), d.g.fam, na.action = na.exclude))
# strong relationship between SLA and days to flowering

p3 = ggplot(d.traits.all, aes(x=BIO12, y=log(SLA_mean), color=region)) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=12), legend.position="none") +
  geom_point(shape=16) +     
  geom_point(shape=1) +    # Use hollow circles
  scale_colour_hue(l=50) + labs(x="mean annual precip. (mm)", y="log(SLA mm2/mg)") +
  geom_smooth(method=lm)   # Add linear regression line 
#  (by default includes 95% confidence region)

p4 = ggplot(d.traits.all, aes(x=BIO12, y=days_to_flower_mean, color=region)) + guides(color=FALSE) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=12), legend.position="none") +
  geom_point(shape=16) +   
  scale_colour_hue(l=50) + labs(x="mean annual precip. (mm)", y="days to flowering") +
  geom_smooth(method=lm, se=T)   # Add linear regression line 
#  (by default includes 95% confidence region)

p5 = ggplot(d.traits.all, aes(x=soil_pc_CV, y=area_plast_nut, color=region)) + guides(color=FALSE) +
  theme_bw() + theme(axis.title=element_text(size=16), axis.text=element_text(size=12), legend.position="none") +
  geom_point(shape=16) +   
  scale_colour_hue(l=50) + labs(x="CV soil nutrients", y="leaf area plasticity to nutrients") +
  geom_smooth(method=lm, se=T)   # Add linear regression line 
#  (by default includes 95% confidence region)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p1, p3, p2, p4, cols=2)

# In general, what we see is an integration of rapid life cycle that is linked to water stress, and also potentially soil nutrients. The set of traits associated with low precipitation (thus short growing season) includes high SLA, fast growth, early flowering, and reduced plasticity in flowering (flower early even at low water supply). 


################### 
# Trait associations

# Make Correlation Figure (Fig 5)

# create matrix of traits for regions separately, and combined
traits = names(d.traits.all)[3:15]
traits_all = d.traits.all[,traits]
traits_CA = filter(traits_all, region =="California")
traits_Chile = filter(traits_all, region =="Chile")
# rename columns for plotting
names(traits_all) = names(traits_CA) = names(traits_Chile) = c("region", "Days to flower", "Stem length", "SLA", "Leaf area", "Days to flower H20 plast", "Stem length H2O plast", "SLA H2O plast", "Leaf Area H2O plast", "Days to flower nut plast", "Stem length nut plast", "SLA nut plast", "Leaf Area nut plast")

cor.all = cor(traits_all[,2:ncol(traits_all)], use="complete.obs")
cor.CA = cor(traits_CA[,2:ncol(traits_CA)], use="complete.obs")
cor.Chile = cor(traits_Chile[,2:ncol(traits_Chile)], use="complete.obs")

res.all <- cor.mtest(traits_all[,2:ncol(traits_all)], 0.05)
res.CA <- cor.mtest(traits_CA[,2:ncol(traits_CA)], 0.05)  
res.Chile <- cor.mtest(traits_Chile[,2:ncol(traits_Chile)], 0.05)

# For plotting, set the diagonal values to 0's -- we want to expand the color scale to be more sensitive, rather than having it necessarily reach all the way to 1.0. Also, the large diagonal blue circles are distracting. 
diag(cor.all) = 0
diag(cor.CA) = 0
diag(cor.Chile) = 0

# Make the plots of correlation matrices & their significance values
par(mar=rep(4, 4))
#pdf("./Plots/corrplot_all.pdf")
corrplot(cor.all, method="circle", order="original", p.mat=res.all[[1]], sig.level=0.05, type="upper", tl.pos="td", is.corr = F, diag=FALSE, cl.lim=c(-0.85, 0.85))
#dev.off()
corrplot(cor.CA, method="circle", order="original", p.mat=res.CA[[1]], sig.level=0.05, type="upper", tl.pos="td", , is.corr = F, diag=FALSE, cl.lim=c(-0.85, 0.85))
corrplot(cor.Chile, method="circle", order="original", p.mat=res.Chile[[1]], sig.level=0.05, type="upper", tl.pos="td", is.corr = F, diag=FALSE, cl.lim=c(-0.85, 0.85))

###############################################################
## END OF ANALYSIS FOR PAPER
###############################################################



