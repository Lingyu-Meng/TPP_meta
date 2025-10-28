# Installing required packages
# install.packages("metafor")
# install.packages("clubSandwich")
# install.packages("MAd")
# install.packages("dplyr")
# install.packages("weightr")
# install.packages("puniform")
# install.packages("metapower")

# Loading required packages

library("metafor")
library("clubSandwich")
library("MAd")
library("dplyr")
library("weightr")
library("puniform")
library("metapower")


# reading the data
MA_dat <- read.csv("Raw.csv", stringsAsFactors= TRUE, encoding="UTF-8")


# data overview

table(MA_dat$IDpaper)
table(table(MA_dat$IDpaper))
length(unique(MA_dat$IDpaper))

table(MA_dat$IDstudy)
table(table(MA_dat$IDstudy))
length(unique(MA_dat$IDstudy))

table(MA_dat$IDsubsample)
table(table(MA_dat$IDsubsample))
length(unique(MA_dat$IDsubsample))

table(MA_dat$IDeffect)
table(table(MA_dat$IDeffect))
length(unique(MA_dat$IDeffect))

table(MA_dat$IDsubsample, MA_dat$DV)
margin.table(table(MA_dat$IDsubsample, MA_dat$DV), 2)



# saving ID variables as factors
MA_dat$IDpaper <- as.factor(MA_dat$IDpaper)
MA_dat$IDstudy <- as.factor(MA_dat$IDstudy)
MA_dat$IDsample <- as.factor(MA_dat$IDsample)
MA_dat$IDsubsample <- as.factor(MA_dat$IDsubsample)
MA_dat$IDeffect <- as.factor(MA_dat$IDeffect)



# creating the covariance matrix for effect sizes from the same subsample
# correlation between same-sample ES is set to .6

V <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.6)


# overall effect size
# fitting the five-level random effects model 

es <- rma.mv(ES_g, V, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es


# checking likelihood profile plots
# par(mar=c(1,1,1,1)) # if plot not showing in RStudio
profile(es)
dev.off()
# par(mar=c(5.1,4.1,4.1,2.1)) # returning the margin values to default (if previously changed)


# obtaining robust variance estimates
coef_test(es, vcov = "CR2")
conf_int(es, vcov = "CR2")

# calculating the overall amount of heterogeneity
sum(es$sigma2)

# obtaining I2 value
W <- diag(1/MA_dat$SV)
X <- model.matrix(es)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es$sigma2) / (sum(es$sigma2) + (es$k-es$p)/sum(diag(P))) # overall I2
100 * es$sigma2 / (sum(es$sigma2) + (es$k-es$p)/sum(diag(P)))  # I2 per levels



# forest plot
forest(es, xlab = "Effect Size", slab = paste (MA_dat$table_label), 
       psize = 1.5, xlim = c(-6, 7), order = order(MA_dat$ES_g), cex = .75, cex.lab = .8)

#saving forest plot to a file
tiff("Forest.tiff", height = 12, width = 8, units = 'in', res = 300)
forest(es, xlab = "Effect Size", slab = paste (MA_dat$table_label), 
       psize = 1.5, xlim = c(-6, 7), order = order(MA_dat$ES_g), cex = .75, cex.lab = .8)
dev.off()


# post-hoc power analysis 

# calculating the number of participants per group
# for between-subject designs N per group = Ntotal / 2
# for within-subject designs N per group = Ntotal

MA_dat$Ngroup <- ifelse(MA_dat$SS == "BS", MA_dat$Ntotal / 2, MA_dat$Ntotal)
MA_dat$Ngroup
MA_dat$Ntotal
mean(MA_dat$Ngroup)

power <- mpower(effect_size = .22, study_size = 82.29, k = 69, i2 = .85, es_type = "d") 
power
plot_mpower(power)



# assessing publication bias / small-study effect

# funnel plot
metafor::funnel(es)

# contour-enhanced funnel plot
metafor::funnel(es, refline=0, level=c(95, 99), 
       shade=c("white", "grey"))

#saving funnel plot to a file
tiff("Funnel.tiff", height = 4, width = 5, units = 'in', pointsize = 7, res = 300)
metafor::funnel(es, refline=0, level=c(95, 99), 
       shade=c("white", "grey"))
dev.off()

# Egger's regression
egger <- rma.mv(ES_g, V, mods = ~ sqrt(SV), random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
egger
coef_test(egger,vcov = "CR2")


# sensitivity models after aggregating and sampling effect sizes

# creating a dataset with effect sizes aggregated per study
MA_dat_agg <- agg(IDpaper, ES_g, SV, n.1=NULL, n.2=NULL, method = "BHHR", cor = .50,  mod=NULL, MA_dat) 

# creating a dataset with a random sampled effect size per study
set.seed(20210204)
MA_dat_sample <- sample_n(group_by(MA_dat, IDpaper), size = 1) 

# Vevea and Hedges 
es.weight_agg <- weightfunct(MA_dat_agg$es, MA_dat_agg$var)
es.weight_agg
es.weight_sample <- weightfunct(MA_dat_sample$ES, MA_dat_sample$SV)
es.weight_sample

# p-uniform*
puni_star(yi = MA_dat_agg$es, vi = MA_dat_agg$var, alpha = 0.05, side = "right")
puni_star(yi = MA_dat_sample$ES, vi = MA_dat_sample$SV, alpha = 0.05, side = "right")


# outliers and influential cases

rstudent(es) # studentised residuals - indicating outliers
# one ES with residual z > 3 can be considered an outlier (Chiang & Wu)
cook <- cooks.distance(es) # indicating influential data points
plot(cook, type="o")
cook
round(cook, digits = 2)
mean(MA_dat$Ntotal)
sort(MA_dat$Ntotal)
# the same study (Chiang & Wu) is an influential case


# overall effect size excluding the influential/outlier effect size (study 09 Chiang & Wu)

# excluding the outlier effect size
MA_dat_s1 <- subset(MA_dat, IDpaper != "9",
                    select=IDpaper:ES_quality) 
# correcting the number of levels in ID factor variables
MA_dat_s1$IDpaper <- as.factor(as.numeric(MA_dat_s1$IDpaper))
MA_dat_s1$IDstudy <- as.factor(as.numeric(MA_dat_s1$IDstudy))
MA_dat_s1$IDsample <- as.factor(as.numeric(MA_dat_s1$IDsample))
MA_dat_s1$IDsubsample <- as.factor(as.numeric(MA_dat_s1$IDsubsample))
MA_dat_s1$IDeffect <- as.factor(as.numeric(MA_dat_s1$IDeffect))

V_s1 <- impute_covariance_matrix(MA_dat_s1$SV, cluster=MA_dat_s1$IDsubsample, r = 0.6)
es_s1 <- rma.mv(ES_g, V_s1, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat_s1)
es_s1
coef_test(es_s1, vcov = "CR2")
conf_int(es_s1, vcov = "CR2")

W <- diag(1/MA_dat_s1$SV)
X <- model.matrix(es_s1)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es_s1$sigma2) / (sum(es_s1$sigma2) + (es_s1$k-es_s1$p)/sum(diag(P)))
100 * es_s1$sigma2 / (sum(es_s1$sigma2) + (es_s1$k-es$p)/sum(diag(P)))




# Single moderator analyses
# models with random effects on the ES level (IDeffect) 
# to obtain correct p-values for moderators intercepts are included in the model
# to obtain effect size values for each level of a categorical moderator, intercepts are excluded (mods = ~ factor-1)


# Participant age
# for moderator analyses regarding age studies Chiang & Wu and Do & Telzer were excluded due to very large age ranges they reported

# excluding the two effect sizes
MA_dat_s2 <- subset(MA_dat, IDpaper != "9" & IDpaper != "18",
                           select=IDpaper:ES_quality)
# correcting the number of levels in ID factor variables
MA_dat_s2$IDpaper <- as.factor(as.numeric(MA_dat_s2$IDpaper))
MA_dat_s2$IDstudy <- as.factor(as.numeric(MA_dat_s2$IDstudy))
MA_dat_s2$IDsample <- as.factor(as.numeric(MA_dat_s2$IDsample))
MA_dat_s2$IDsubsample <- as.factor(as.numeric(MA_dat_s2$IDsubsample))
MA_dat_s2$IDeffect <- as.factor(as.numeric(MA_dat_s2$IDeffect))

V_s2 <- impute_covariance_matrix(MA_dat_s2$SV, cluster=MA_dat_s2$IDsubsample, r = 0.6)

# plotting the age-ES relationship
plot(MA_dat_s2$NewAge, MA_dat_s2$ES_g, pch=19, cex=.14/sqrt(MA_dat_s2$SV),
     xlab="Age", ylab="Effect Size", main="Scatter plot")

#saving scatter plot to a file
tiff("Scatter.tiff", height = 4, width = 5, units = 'in', pointsize = 7, res = 300)
plot(MA_dat_s2$NewAge, MA_dat_s2$ES_g, pch=19, cex=.14/sqrt(MA_dat_s2$SV),
     xlab="Age", ylab="Effect Size", main="Scatter plot")
dev.off()

# linear model
es.NewAge <- rma.mv(ES_g, V_s2, 
                    mods = ~ NewAge,
                    random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat_s2)
es.NewAge 
coef_test(es.NewAge, vcov = "CR2")
conf_int(es.NewAge, vcov = "CR2")
W <- diag(1/MA_dat_s2$SV)
X <- model.matrix(es.NewAge)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.NewAge$sigma2) / (sum(es.NewAge$sigma2) + (es.NewAge$k-es.NewAge$p)/sum(diag(P)))


# Game type and interdependence

# for the moderator analysis regarding game type Do & Telzer study was excluded 
# in order for all game type moderator levels to be represented by more than 1 study

# excluding the effect size
MA_dat_s3 <- subset(MA_dat, IDpaper != "18",
                    select=IDpaper:ES_quality)
# correcting the number of levels in ID factor variables
MA_dat_s3$IDpaper <- as.factor(as.numeric(MA_dat_s3$IDpaper))
MA_dat_s3$IDstudy <- as.factor(as.numeric(MA_dat_s3$IDstudy))
MA_dat_s3$IDsample <- as.factor(as.numeric(MA_dat_s3$IDsample))
MA_dat_s3$IDsubsample <- as.factor(as.numeric(MA_dat_s3$IDsubsample))
MA_dat_s3$IDeffect <- as.factor(as.numeric(MA_dat_s3$IDeffect))

V_s3 <- impute_covariance_matrix(MA_dat_s3$SV, cluster=MA_dat_s3$IDsubsample, r = 0.6)


es.DV <- rma.mv(ES_g, V_s3, 
                mods = ~ DV,
                random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat_s3)
es.DV 
coef_test(es.DV, vcov = "CR2")
conf_int(es.DV, vcov = "CR2")
W <- diag(1/MA_dat_s3$SV)
X <- model.matrix(es.DV)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.DV$sigma2) / (sum(es.DV$sigma2) + (es.DV$k-es.DV$p)/sum(diag(P)))



es.InterDep <- rma.mv(ES_g, V, 
                mods = ~ InterDep,
                random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es.InterDep 
coef_test(es.InterDep, vcov = "CR2")
conf_int(es.InterDep, vcov = "CR2")
W <- diag(1/MA_dat$SV)
X <- model.matrix(es.InterDep)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.InterDep$sigma2) / (sum(es.InterDep$sigma2) + (es.InterDep$k-es.InterDep$p)/sum(diag(P)))

# a stricter test of the interdependence moderator
# analysis on experimentally created groups only

# selecting experimentally created groups only
MA_dat_s4 <- subset(MA_dat, IV == "E",
                    select=IDpaper:ES_quality)
# correcting the number of levels in ID factor variables
MA_dat_s4$IDpaper <- as.factor(as.numeric(MA_dat_s4$IDpaper))
MA_dat_s4$IDstudy <- as.factor(as.numeric(MA_dat_s4$IDstudy))
MA_dat_s4$IDsample <- as.factor(as.numeric(MA_dat_s4$IDsample))
MA_dat_s4$IDsubsample <- as.factor(as.numeric(MA_dat_s4$IDsubsample))
MA_dat_s4$IDeffect <- as.factor(as.numeric(MA_dat_s4$IDeffect))

V_s4 <- impute_covariance_matrix(MA_dat_s4$SV, cluster=MA_dat_s4$IDsubsample, r = 0.6)

es.InterDep_s4 <- rma.mv(ES_g, V_s4, 
                      mods = ~ InterDep,
                      random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat_s4)
es.InterDep_s4 
coef_test(es.InterDep_s4, vcov = "CR2")
conf_int(es.InterDep_s4, vcov = "CR2")
W <- diag(1/MA_dat_s4$SV)
X <- model.matrix(es.InterDep_s4)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.InterDep_s4$sigma2) / (sum(es.InterDep_s4$sigma2) + (es.InterDep_s4$k-es.InterDep_s4$p)/sum(diag(P)))


# Membership manipulation, group type, currency type and sex

es.SS <- rma.mv(ES_g, V, 
             mods = ~ SS,
             random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es.SS 
coef_test(es.SS, vcov = "CR2")
conf_int(es.SS, vcov = "CR2")
W <- diag(1/MA_dat$SV)
X <- model.matrix(es.SS)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.SS$sigma2) / (sum(es.SS$sigma2) + (es.SS$k-es.SS$p)/sum(diag(P)))


es.IV <- rma.mv(ES_g, V, 
                mods = ~ IV,
                random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es.IV 
coef_test(es.IV, vcov = "CR2")
conf_int(es.IV, vcov = "CR2")
W <- diag(1/MA_dat$SV)
X <- model.matrix(es.IV)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.IV$sigma2) / (sum(es.IV$sigma2) + (es.IV$k-es.IV$p)/sum(diag(P)))


es.CurrencyType <- rma.mv(ES_g, V, 
                mods = ~ CurrencyType,
                random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es.CurrencyType 
coef_test(es.CurrencyType, vcov = "CR2")
conf_int(es.CurrencyType, vcov = "CR2")
W <- diag(1/MA_dat$SV)
X <- model.matrix(es.CurrencyType)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.CurrencyType$sigma2) / (sum(es.CurrencyType$sigma2) + (es.CurrencyType$k-es.CurrencyType$p)/sum(diag(P)))



which(is.na(MA_dat$sex))
MA_dat$IDeffect[which(is.na(MA_dat$sex))]
MA_dat$IDstudy[which(is.na(MA_dat$sex))]

# for the moderator analysis regarding sex Bauer et al. and List et al. studies were excluded
# due to missing values for proportion of males

# excluding the two effect sizes
MA_dat_s5 <- subset(MA_dat, IDpaper != "8" & IDpaper != "13",
                    select=IDpaper:ES_quality)
# correcting the number of levels in ID factor variables
MA_dat_s5$IDpaper <- as.factor(as.numeric(MA_dat_s5$IDpaper))
MA_dat_s5$IDstudy <- as.factor(as.numeric(MA_dat_s5$IDstudy))
MA_dat_s5$IDsample <- as.factor(as.numeric(MA_dat_s5$IDsample))
MA_dat_s5$IDsubsample <- as.factor(as.numeric(MA_dat_s5$IDsubsample))
MA_dat_s5$IDeffect <- as.factor(as.numeric(MA_dat_s5$IDeffect))

V_s5 <- impute_covariance_matrix(MA_dat_s5$SV, cluster=MA_dat_s5$IDsubsample, r = 0.6)


es.sex <- rma.mv(ES_g, V_s5, 
                          mods = ~ sex,
                          random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat_s5)
es.sex 
coef_test(es.sex, vcov = "CR2")
conf_int(es.sex, vcov = "CR2")
W <- diag(1/MA_dat_s5$SV)
X <- model.matrix(es.sex)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.sex$sigma2) / (sum(es.sex$sigma2) + (es.sex$k-es.sex$p)/sum(diag(P)))


# Country

table(MA_dat$country)

es.country <- rma.mv(ES_g, V, 
                    mods = ~ country-1,
                    random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es.country 
coef_test(es.country, vcov = "CR2")
conf_int(es.country, vcov = "CR2")
W <- diag(1/MA_dat$SV)
X <- model.matrix(es.country)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(es.country$sigma2) / (sum(es.country$sigma2) + (es.country$k-es.country$p)/sum(diag(P)))






################################################################################################

# sensitivity analysis for the value of correlation of effect sizes obtained on the same group of participants


V00 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0)
V10 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.1)
V20 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.2)
V30 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.3)
V40 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.4)
V50 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.5)
V60 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.6)
V70 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.7)
V80 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.8)
V90 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 0.9)
V100 <- impute_covariance_matrix(MA_dat$SV, cluster=MA_dat$IDsubsample, r = 1)


es00 <- rma.mv(ES_g, V00, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es10 <- rma.mv(ES_g, V10, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es20 <- rma.mv(ES_g, V20, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es30 <- rma.mv(ES_g, V30, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es40 <- rma.mv(ES_g, V40, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es50 <- rma.mv(ES_g, V50, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es60 <- rma.mv(ES_g, V60, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es70 <- rma.mv(ES_g, V70, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es80 <- rma.mv(ES_g, V80, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es90 <- rma.mv(ES_g, V90, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)
es100 <- rma.mv(ES_g, V100, random = ~ 1 | IDpaper / IDstudy / IDsubsample / IDeffect, data=MA_dat)

es00 # es = .2116
es10 # es = .2133
es20 # es = .2151
es30 # es = .2170
es40 # es = .2193
es50 # es = .2221
es60 # es = .2229
es70 # es = .2119
es80 # es = .2204
es90 # es = .2186
es100 # es = .2116
