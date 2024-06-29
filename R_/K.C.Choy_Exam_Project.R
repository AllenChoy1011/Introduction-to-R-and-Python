# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Examination project in R 
# Author: Kin Cheung Choy
# Email: allenreve@icloud.com
# Submission date:04012022
# Version: 1
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  


# Data Management ---------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: data import, variable assignment, dataset reorganisation (merge + long format),
## Load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
library(GGally)
library(SciViews)

## set work directory 
getwd()
setwd('~/Desktop/R')

## data import
data_pk <- read.table(file = "BPI889_PK_55.csv", header = T, sep = ',', na.strings = '.')
data_snp <- read.table(file = "BPI889_SNP_55.txt", header = T, row.names = NULL)

## Visually inspect imported data.frames
str(data_pk)
str(data_snp)

## Variable assignment
### Set variables names
names(data_pk)[1] <- 'ID'
names(data_snp)[1] <- 'ID'

## Dataset reorganisation
### Merge datasets
data_merge <- merge(data_pk, data_snp, by = "ID")

### Transform dataset from wide to long format
data_tidy <- gather(data_merge, Time, Con, 2:16)
data_all <- data_tidy %>% select(ID, Time, Con, T134A, A443G, G769C, G955C, A990C)

### Convert variable into numberic
data_all$Time <- as.numeric(gsub('Time.|.h', '', data_all$Time))

### Convert variable into factor
data_all$ID <- as.factor(data_all$ID)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  

# Variable calculations ---------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: calculation of body size measurement, categorization of body size measurement, 
# PK variable calculation

## Body size measurement(TBW)
MTBW <- 2.447 - (0.09156 * data_pk$Age..yrs.) + (0.1074 * data_pk$Height..cm.) + (0.3362 * data_pk$Weight..kg.)
FTBW <- -2.097 + (0.1069 * data_pk$Height..cm.) + (0.2466 * data_pk$Weight..kg.)
data_all$TBW <- ifelse(data_pk$Sex == "M", MTBW, FTBW)
data_all$TBW <- round(data_all$TBW, digits = 1)

## Categorization of body size measurement
data_all$CTBW <- ifelse(data_all$TBW > 40, "H", "L")

## Cmax Calculation
data_all$Cmax <- NA
data_all$Con <- as.numeric(data_all$Con)
data_all <- data_all %>% group_by(ID) %>% mutate(Cmax = max(Con, na.rm = T))

## Vd Calculation
### Add new variable lnC
data_all$lnC <- ln(data_all$Con)
data_all$lnC[is.infinite(data_all$lnC)] <- NA 

### Add linear regression 
data_filter <- filter(data_all, Time > 5 & Time < 24)
model_interact <- lm(formula = lnC ~ ID * Time + 0, data_filter) # x + 0 to exclude the intercept 
summary(model_interact)

### lnC0 calculation
lnC0_var <- model_interact$coefficients[1:100] %>%
  as.data.frame()
lnC0_var$ID <- rownames(lnC0_var)
lnC0_var$ID <- sub("..", "", lnC0_var$ID)
colnames(lnC0_var)[1] <- "lnC0"

### k constant calculation
k_var <- model_interact$coefficients[101:200] %>% 
  as.data.frame()
colnames(k_var)[1] <- "k"
k_var$`k`[2:100] <- k_var$`k`[2:100] + (-0.74)
k_var$ID <- rownames(k_var)
k_var$ID <- gsub(".{5}$", "", k_var$ID )
k_var$ID <- sub("..", "", k_var$ID )
k_var$ID[1] <- "pat1"

### Merge k and lnC0 
data_all <- merge(data_all, lnC0_var, by="ID")
data_all <- merge(data_all, k_var, by="ID")

### Vd calculation
data_all$Vd <- 200/exp(data_all$lnC0)    

### CL calculation
data_all$CL <- (-data_all$k) * data_all$Vd 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  

# Data Exploration --------------------------------------------------------
summary(data_all)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: numerical summary of PK variables, graphical assessment of 1) PK profiles,
# 2) PK variable correlations, 3)PK variable-SNP correlations, 
# 4) PK variable-body size measurement correlation with linear regression
###individual concentrations of BPI889 versus time (spaghetti plot)
p1 <- ggplot(data = data_all, aes(x = Time, y = Con, group = ID, shape = ID, col = ID)) + 
  geom_line(size = 0.5) +
  xlab("Time") + ylab("Concentration(mg/L)") + 
  ggtitle("Individual concentrations of BPI889 versus time ") +
  theme_bw()
p1
###correlations between Cmax,Vd and CL (scatter plot)
p2 <- pairs(~ Cmax + Vd + CL, data = data_all,
            upper.panel = panel.smooth,
            pch = 19,
            col = 'blue',
            bg = 'blue')
            
p2 
###correlation between Vd and SNPs 
p3 <- list()
for (i in colnames(data_all[,4:8])) {
  p3[[i]] <- ggplot(data_all, aes(x = factor(data_all[[i]]), 
                                  y = Vd, 
                                  fill = factor(data_all[[i]]))) +
    geom_boxplot(outlier.colour = "grey", outlier.shape = 14,
                 outlier.size = 1, notch = T,
                 show.legend = F) +
    labs(x = i) + 
    labs(title = "Box plot of Vd by SNPs")
}

ggarrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],p3[[5]])
###correlation between CL and SNPs
p4 <- list()
for (i in colnames(data_all[,4:8])) {
  p4[[i]] <- ggplot(data_all, aes(x = factor(data_all[[i]]), 
                                  y = CL, 
                                  fill = factor(data_all[[i]]))) +
    geom_boxplot(outlier.colour = "grey", outlier.shape = 14,
                 outlier.size = 1, notch = T,
                 show.legend = F) +
    labs(x = i) + 
    labs(title = "Box plot of CL by SNPs")
}

ggarrange(p4[[1]],p4[[2]],p4[[3]],p4[[4]],p[[5]])

###correlations between Vd and TBW to assess a relationship and add a linear regression (scatter with linear regression)
p5 <- ggpairs(data = data_all, 
              columns = c(9,15),
              lower = list(continuous = "smooth"))
p5
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  



# Statistical testing -----------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
# Must include: ANOVA of PK variables for SNPs, t-test of PK variable for body size measurement groups 
###ANOVA of CL for the five SNPs
for (i in colnames(data_all[,4:8])) {
  CL_SNP.aov <- lm(CL ~ data_all[[i]] + 0, data = data_all)
}
summary(CL_SNP.aov)
###ANOVA of Cmax for the five SNPs
for (i in colnames(data_all[,4:8])) {
  Cmax_SNP.aov <- lm(Cmax ~ data_all[[i]] + 0, data = data_all)
}
summary(Cmax_SNP.aov)
###t-test of Vd for the two categorical groups of TBW
t.test(data = data_all, Vd ~ CTBW)
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---  
