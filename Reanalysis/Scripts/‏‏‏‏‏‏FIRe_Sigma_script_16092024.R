####### Eilat photo-physiology
# Loading data (Tamar @ home): 
setwd("C:/Users/Lenovo/Documents/STUDY/Data_Analayzing/Data_CSVs")
FIRe_ave<- read_csv('FIRe_ave.csv')
str(FIRe_ave)

# Packages to call from library when opening R:
# (to add new ones to this list)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(car)
library(tidyr)
library(emmeans)
library(reshape2)
library(scales)
library(ARTool)
library(multcomp)
library(rcompanion)
library(psych)

# If we want to instruct R to avoid using scientific notation and instead 
# display numbers in their full numeric form we use: options(scipen = 999)
options(scipen = 0) # To display numbers in scientific notation.

##### Setting the DEPTH, SITE & Species as factors: 
FIRe_ave$DEPTH <- as_factor(FIRe_ave$DEPTH) 
FIRe_ave$SITE <- as_factor(FIRe_ave$SITE) 
FIRe_ave$Species <- as_factor(FIRe_ave$Species) 

#### Removing all rows with NA's
FIRe_ave <- FIRe_ave %>% na.omit() 
str(FIRe_ave)

################ Normality tests for the Sigma data: 

#### Kolmogorov-Smirnov test:

# K-S test for SITE groups:
ks.test(FIRe_ave[FIRe_ave$SITE=="IUI_S",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$SITE=="IUI_N",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$SITE=="JG",]$Sigma, "pnorm")

# K-S test for DEPTH groups:
ks.test(FIRe_ave[FIRe_ave$DEPTH=="45",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$DEPTH=="35",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$DEPTH=="25",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$DEPTH=="15",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$DEPTH=="5",]$Sigma, "pnorm")

# K-S test for Species groups:
ks.test(FIRe_ave[FIRe_ave$Species=="Stylophora",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$Species=="Porites",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$Species=="Pocillopora",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$Species=="Leptoseris",]$Sigma, "pnorm")
ks.test(FIRe_ave[FIRe_ave$Species=="Alveopora",]$Sigma, "pnorm")


#################### Q-1 (by coral): For each coral what is the optimal depth?

# Subset the data for the Species (to 5 different DF's):
FIRe_ave_Porites <- subset(FIRe_ave, Species == "Porites")
FIRe_ave_Stylophora <- subset(FIRe_ave, Species == "Stylophora")
FIRe_ave_Pocillopora <- subset(FIRe_ave, Species == "Pocillopora")
FIRe_ave_Leptoseris <- subset(FIRe_ave, Species == "Leptoseris")
FIRe_ave_Alveopora <- subset(FIRe_ave, Species == "Alveopora")


############### ART-Anova test for the Sigma parameter for the 5 corals:
# Required packages: install.packages("ARTool"), install.packages("emmeans"), install.packages("multcomp")
# and: install.packages("rcompanion"), install.packages("ggplot2"), install.packages("psych")

################### STYLOPHORA:
### Aligned ranks anova:
model_Stylo_Sigma = art(Sigma ~ DEPTH+(1|SITE) , data = FIRe_ave_Stylophora) # +(1|SITE) took off the site factor
### Check the success of the procedure:
model_Stylo_Sigma

### Conduct ANOVA:
anova(model_Stylo_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Stylo_DEPTH_posthoc = art.con(model_Stylo_Sigma, "DEPTH", adjust="tukey")
Stylo_DEPTH_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Stylo_DEPTH = as.data.frame(Stylo_DEPTH_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Stylo_DEPTH)


################### PORITES:
### Aligned ranks anova:
model_Porites_Sigma = art(Sigma ~ DEPTH+(1|SITE) , data = FIRe_ave_Porites)
### Check the success of the procedure:
model_Porites_Sigma

### Conduct ANOVA:
anova(model_Porites_Sigma)  ## Non significant!


################### POCILLOPORA:
### Aligned ranks anova:
model_Pocill_Sigma = art(Sigma ~ DEPTH+(1|SITE) , data = FIRe_ave_Pocillopora)
### Check the success of the procedure:
model_Pocill_Sigma

### Conduct ANOVA:
anova(model_Pocill_Sigma)  ## Non significant!


################### LEPTOSERIS:
### Aligned ranks anova:
model_Lepto_Sigma = art(Sigma ~ DEPTH+(1|SITE) , data = FIRe_ave_Leptoseris)
### Check the success of the procedure:
model_Lepto_Sigma

### Conduct ANOVA:
anova(model_Lepto_Sigma) ## Non significant! 


################### ALVEOPORA:
### Aligned ranks anova:
model_Alveo_Sigma = art(Sigma ~ DEPTH , data = FIRe_ave_Alveopora)
### Check the success of the procedure:
model_Alveo_Sigma

### Conduct ANOVA:
anova(model_Alveo_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Alveo_DEPTH_posthoc = art.con(model_Alveo_Sigma, "DEPTH", adjust="tukey")
Alveo_DEPTH_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Alveo_DEPTH = as.data.frame(Alveo_DEPTH_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Alveo_DEPTH)


############# Plotting Sigma: 

# Calculates the 25th percentile=מאון (1st quartile) and the 75th percentile (3rd quartile) 
# of the Sigma column in the FIRe_ave DF:
quartiles <- quantile(FIRe_ave$Sigma, probs = c(.25,.75))
# Calculates the interquartile range (IQR) of the Sigma column 
# which is the difference between the 3rd and 1st quartiles:
IQR_cor <- IQR(FIRe_ave$Sigma)

# To define the lower & the upper bound for detecting outliers:
lower <- quartiles[1] - 1.5*IQR_cor # calculated as the first quartile minus 1.5 times the IQR
upper <- quartiles[2] + 1.5*IQR_cor # calculated as the first quartile plus 1.5 times the IQR

# To create new DF with the rows that values in it are within bounds above (only), without outliers:
Sigma_no_outlier <- subset(FIRe_ave, FIRe_ave$Sigma >lower & FIRe_ave$Sigma< upper)

Sigma <- ggplot(data = Sigma_no_outlier,aes(x= DEPTH, y=Sigma, fill= Species))+
  geom_boxplot() + theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(x="Depth in meters", y="Sigma (A^2)") + # axes labels 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(colour="Black", size=16),
        axis.title.y = element_text(colour="Black", size=16),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        # strip.background = element_blank(),  # Optional: remove background of facet labels
        strip.text = element_text(size = 13) # Customize facet label text size
  ) +
  # axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev) + # Depth from shallow to deep
  coord_flip() + # Flips the coordinates, making the x-axis vertical and the y-axis horizontal.
  facet_grid(Species~., scales="free", space="free")

# Print the plot
print(Sigma)   


################### Q-2 (by depth): For each depth level, which is the best preforming coral?

# Subset the data for the Species (to 5 different DF's):
FIRe_ave_5 <- subset(FIRe_ave, DEPTH == "5")
FIRe_ave_15 <- subset(FIRe_ave, DEPTH == "15")
FIRe_ave_25 <- subset(FIRe_ave, DEPTH == "25")
FIRe_ave_35 <- subset(FIRe_ave, DEPTH == "35")
FIRe_ave_45 <- subset(FIRe_ave, DEPTH == "45")

############### ART-Anova test for the Sigma parameter for the 5 depth-points:

################### 5m:
### Aligned ranks anova:
model_5_Sigma = art(Sigma ~ Species+(1|SITE) , data = FIRe_ave_5) # +(1|SITE) took off the site factor
### Check the success of the procedure:
model_5_Sigma

### Conduct ANOVA:
anova(model_5_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Species_5_posthoc = art.con(model_5_Sigma, "Species", adjust="tukey")
Species_5_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Species_5 = as.data.frame(Species_5_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Species_5)


################### 15m:
### Aligned ranks anova:
model_15_Sigma = art(Sigma ~ Species+(1|SITE) , data = FIRe_ave_15) # +(1|SITE) took off the site factor
### Check the success of the procedure:
model_15_Sigma

### Conduct ANOVA:
anova(model_15_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Species_15_posthoc = art.con(model_15_Sigma, "Species", adjust="tukey")
Species_15_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Species_15 = as.data.frame(Species_15_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Species_15)


################### 25m:
### Aligned ranks anova:
model_25_Sigma = art(Sigma ~ Species, data = FIRe_ave_25) 
### Check the success of the procedure:
model_25_Sigma

### Conduct ANOVA:
anova(model_25_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Species_25_posthoc = art.con(model_25_Sigma, "Species", adjust="tukey")
Species_25_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Species_25 = as.data.frame(Species_25_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Species_25)


################### 35m:
### Aligned ranks anova:
model_35_Sigma = art(Sigma ~ Species, data = FIRe_ave_35) 
### Check the success of the procedure:
model_35_Sigma

### Conduct ANOVA:
anova(model_35_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Species_35_posthoc = art.con(model_35_Sigma, "Species", adjust="tukey")
Species_35_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Species_35 = as.data.frame(Species_35_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Species_35)


################### 45m:
### Aligned ranks anova:
model_45_Sigma = art(Sigma ~ Species, data = FIRe_ave_45) 
### Check the success of the procedure:
model_45_Sigma

### Conduct ANOVA:
anova(model_45_Sigma)

# Posthoc pairwise comparison for the DEPTH:
Species_45_posthoc = art.con(model_45_Sigma, "Species", adjust="tukey")
Species_45_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_Species_45 = as.data.frame(Species_45_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_Species_45)


############# Plotting Sigma by facet depth:

# Calculates the 25th percentile=מאון (1st quartile) and the 75th percentile (3rd quartile) 
# of the Sigma column in the FIRe_ave DF:
quartiles <- quantile(FIRe_ave$Sigma, probs = c(.25,.75))
# Calculates the interquartile range (IQR) of the Sigma column 
# which is the difference between the 3rd and 1st quartiles:
IQR_cor <- IQR(FIRe_ave$Sigma)

# To define the lower & the upper bound for detecting outliers:
lower <- quartiles[1] - 1.5*IQR_cor # calculated as the first quartile minus 1.5 times the IQR
upper <- quartiles[2] + 1.5*IQR_cor # calculated as the first quartile plus 1.5 times the IQR

# To create new DF with the rows that values in it are within bounds above (only), without outliers:
Sigma_no_outlier <- subset(FIRe_ave, FIRe_ave$Sigma >lower & FIRe_ave$Sigma< upper)

Sigma <- ggplot(data = Sigma_no_outlier,aes(x= Species, y=Sigma, fill= Species))+
  geom_boxplot() + theme_classic() + stat_boxplot(geom = 'errorbar', width = 0.75, position = "dodge") +
  labs(x="Corals", y="Sigma (A^2)") + # axes labels 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title.x = element_text(colour="Black", size=16),
        axis.title.y = element_text(colour="Black", size=16),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        # strip.background = element_blank(),  # Optional: remove background of facet labels
        strip.text = element_text(size = 13) # Customize facet label text size
  ) +
  # axis.title.y = element_blank()) +
  scale_x_discrete(limits = rev) + # Depth from shallow to deep
  coord_flip() + # Flips the coordinates, making the x-axis vertical and the y-axis horizontal.
  facet_grid(DEPTH~., scales="free", space="free")

# Print the plot 
print(Sigma)  


