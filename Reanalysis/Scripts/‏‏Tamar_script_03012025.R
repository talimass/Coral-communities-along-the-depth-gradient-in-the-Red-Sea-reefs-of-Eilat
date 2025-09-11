# title: "Distribution and photo-physiology of coral communities across the 
# depth gradient in the Red Sea reefs of Eilat"
# author: "Tamar  Shifroni-Taylor"
# date: "04/04/2024"

  # Packages to call from library when opening R:
  # (to add new ones to this list)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(car)
library(reshape2)
library(scales)
library(ARTool)
library(tidyr)
library(emmeans)
library(multcomp)
library(viridis)
library(rcompanion)
library(psych)
library(devtools)
library(pairwiseAdonis)
### Installed the package 'pairwiseAdonis' with the two command below:
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


# 28/10/24 - I didn't run this one yet (causing trouble...): (Got results from Maya)
# library(LDM) 

# Loading data by Maya @ Uni:
# please note, when you copy a path from windows, the slash is in opposite direction. 
# If not changed, R will not be able to read the path.
# setwd("F:/Tamar_Mass")# set working directory. remove this line from published code
# lc_sf.data<- read.csv('Annotations _MO_04042024_ver1.csv', header=TRUE, stringsAsFactors = TRUE)

# Loading data (Tamar @ home): 
setwd("C:/Users/Lenovo/Documents/STUDY/Data_Analayzing/Data_CSVs")
lc_sf.data<-read.csv('Annotations_MO_25062024_ver1.csv',header=TRUE, stringsAsFactors = TRUE)
str(lc_sf.data)
levels(lc_sf.data$GROUP)

# because we do not really want to tread depth as continuous value, 
# we will change it to be treated as factor:
lc_sf.data$DEPTH=as.factor(lc_sf.data$DEPTH)


################ calculate proportion of live:
# to summarize how many live  and how many Non_live we have per picture:
live=aggregate(Count~SITE+DEPTH+PIC+Live_Non_live, data=lc_sf.data, FUN="sum") 
# to keep only live data:
live=live[live$Live_Non_live=="Live",]

# to calculate the proportion:
live$prop=live$Count/64

# add variable that gives % coverage instead of proportion
live$Cov_percent=live$prop*100
min(live$Cov_percent)
max(live$Cov_percent)
median(live$Cov_percent)

# add variable that combines site and depth
live$site_depth=paste(live$SITE, live$DEPTH, sep="_")


############### Live coverage initial plot (my addition):
ggplot(live, aes(x=DEPTH, y=prop,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + # ggtitle("Live observations proportions in each site by depth") +
 scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Depth (M)") + ylab("Precentage") +
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 

############### Homogeneity of variance Levene's tests:

# install.packages("car")
# library(car)
prop_data<-live$prop # on prop data of all sites & depths (just a trial) 

# For site groups:
leveneTest(prop_data~SITE, data=live)

# For depth groups:
leveneTest(prop_data~DEPTH, data=live)


################ Normality tests: 

#### Kolmogorov-Smirnov test (For larger data, over 50 in a sample):   

# K-S test for SITE groups:
ks.test(live[live$SITE=="IUI_S",]$prop, "pnorm")
ks.test(live[live$SITE=="IUI_N",]$prop, "pnorm")
ks.test(live[live$SITE=="JG",]$prop, "pnorm")

# K-S test for DEPTH groups:
ks.test(live[live$DEPTH=="45",]$prop, "pnorm")
ks.test(live[live$DEPTH=="35",]$prop, "pnorm")
ks.test(live[live$DEPTH=="25",]$prop, "pnorm")
ks.test(live[live$DEPTH=="15",]$prop, "pnorm")
ks.test(live[live$DEPTH=="5",]$prop, "pnorm")

# add variable that gives arcsin transformation instead of proportion:
live$arcsin_prop<-asin(live$prop)

### K-S normality test on the arcsin_prop in the way I did on 'prop':

# K-S test on arcsin prop for SITE groups:
ks.test(live[live$SITE=="IUI_S",]$arcsin_prop, "pnorm")
ks.test(live[live$SITE=="IUI_N",]$arcsin_prop, "pnorm")
ks.test(live[live$SITE=="JG",]$arcsin_prop, "pnorm")

# K-S test on arcsin prop for DEPTH groups:
ks.test(live[live$DEPTH=="45",]$arcsin_prop, "pnorm")
ks.test(live[live$DEPTH=="35",]$arcsin_prop, "pnorm")
ks.test(live[live$DEPTH=="25",]$arcsin_prop, "pnorm")
ks.test(live[live$DEPTH=="15",]$arcsin_prop, "pnorm")
ks.test(live[live$DEPTH=="5",]$arcsin_prop, "pnorm")


############## Density plots for live prop by site:

ggplot(data = live, aes(x = prop, fill = SITE)) +
  geom_density(alpha = 0.7) + # Adjust alpha for transparency if needed
  scale_fill_manual(values = c("JG" = "slategray", "IUI_S"="lightseagreen", "IUI_N" = "navyblue")) + 
  theme_bw() +
  ylim(0,3.8) + # change limit of y axis
  # ggtitle("Density plot of live proportion by site") +
  xlab("Live proportion") + ylab("Density") +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.position = c(1,1), # legend position in the top right of plot area
        legend.justification = c(1,1), # legend position inside the plot area
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

########## Density plots for live prop by depth:

ggplot(data = live, aes(x = prop, fill = DEPTH)) +
  geom_density(alpha = 0.7) + # Adjust alpha for transparency if needed
  scale_fill_manual(values = c("45"="slategray", "35"="lightcoral", "25"="lightseagreen", "15"="maroon", "5"="midnightblue")) + 
  theme_bw() +
  ylim(0,3.8) + # change limit of y axis
  # ggtitle("Density plot of live proportion by depth") +
  xlab("Live proportion") + ylab("Density") +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.position = c(1,1), # legend position in the top right of plot area
        legend.justification = c(1,1), # legend position inside the plot area
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )


################ ART ANOVA test for the 'live' data:

# Required packages: install.packages("ARTool"), install.packages("emmeans"), install.packages("multcomp")
# and: install.packages("rcompanion"), install.packages("ggplot2"), install.packages("psych")

### Aligned ranks anova:
model_live = art(prop ~ DEPTH + SITE + DEPTH:SITE, data = live)
### Check the success of the procedure:
model_live

### Conduct ANOVA:
anova(model_live)

# Posthoc pairwise comparison for SITE:
SITE_posthoc = art.con(model_live, "SITE")
SITE_posthoc
Sum_SITE = as.data.frame(SITE_posthoc) # library(rcompanion)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_SITE)


# Posthoc pairwise comparison for DEPTH:
DEPTH_posthoc = art.con(model_live, "DEPTH")
DEPTH_posthoc
Sum_DEPTH = as.data.frame(DEPTH_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_DEPTH)


# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_posthoc = art.con(model_live, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_DEPTH_SITE = as.data.frame(DEPTH_SITE_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE)


##### Calculating medians of 'live' props for every 'site_depth' group:

# Load dplyr package first: library(dplyr)

# Calculate the median of 'prop' for each 'site_depth' combination
median_site_depth <- live %>%
  group_by(site_depth) %>%
  summarize(median_prop = median(prop, na.rm = TRUE))


##### MAYA_composition ### setwd("F:/Tamar_Mass")

setwd("F:/Tamar_Mass")
data_all=lc_sf.data
write.csv(data_all, "data_all.csv")

Labels_corrected=read.csv("Labels_corrected_for_genus.csv", header=TRUE, stringsAsFactors = TRUE)
str(Labels_corrected)

# Importing corrected labels into data_all DF from Labels_corrected DF:
data_all$LABEL_corrected=Labels_corrected$LABEL1[match(data_all$LABEL,Labels_corrected$LABEL)]
str(data_all)

# Aggregation of observation by label by photo (73 taxons)
data_all_image=aggregate(Count~PIC+LABEL_corrected, data=data_all, FUN="sum")

# Change to % 
data_all_image$coverage_label_cor=data_all_image$Count/64*100

# Adding column present / not present:
data_all_image$prevalence=ifelse(data_all_image$Count>0,1,0)
str(data_all_image)

# Pivot the DF to wide & fill 0 for 0 values using library(reshape2):
data_all_image_wide=as.data.frame(data_all_image[,c(1:2,4)] %>%
pivot_wider(names_from = LABEL_corrected, values_from = coverage_label_cor, values_fill=0))
row.names(data_all_image_wide)=data_all_image_wide$PIC
data_all_image_wide=as.data.frame(t(data_all_image_wide[,-1]))

write.csv(data_all_image_wide, "data_all_image_wide.csv")


############# Create metadata table for images using library(vegan): 
metadata=unique(data_all[,1:3])
write.csv(metadata, "metadata.csv")

# Order the names (photos) in 2 DF's the same by creating the 'ord' vector:
row.names(metadata)=metadata$PIC
ord=match(colnames(data_all_image_wide), row.names(metadata))
metadata=metadata[ord,]

# Hellinger transform the data and create new DF:
data_all_image_wide_hell=decostand(data_all_image_wide, method="total", MARGIN=2)

# Removing the Non_live row to run the test for live only: 
levels(Labels_corrected$LABEL1)
remove_Non_live=as.vector(c("Dead_coral","Other","Rock" ,"Rubble","Sand"))
data_all_image_wide_hell_live=data_all_image_wide_hell[row.names(data_all_image_wide_hell)!="Non_live",]



######## Functional groups plotting using library(ggplot2), library(reshape2):

   ###### Loading data (Tamar @ home desktop): 
setwd("C:/Users/Lenovo/Documents/STUDY/Data_Analayzing/Data_CSVs")
Fun_groups_plot<-read.csv('Fun_groups_plot.csv',header=TRUE, stringsAsFactors = TRUE)
str(Fun_groups_plot) 

# Change depth to factor:
Fun_groups_plot$Depth=as.factor(Fun_groups_plot$Depth)

# Ensure Site factor levels are correct
Fun_groups_plot$Site <- factor(Fun_groups_plot$Site, levels = c("IUI_N", "IUI_S", "JG"))

# Reorder Site levels to "JG" at the top, "IUI_N" in the middle, "IUI_S" at the bottom
Fun_groups_plot$Site <- factor(Fun_groups_plot$Site, levels = c("JG", "IUI_N", "IUI_S"))

# Melt the data frame to long format for ggplot
df_long <- melt(Fun_groups_plot, id.vars = c("Site", "Depth"), variable.name = "Category", value.name = "Percentage")

# Create 100% stacked bar plot:
Groups_plot <- ggplot(df_long, aes(x = Percentage, y = Depth, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) + # to change to stacked plot that not 100% - position = "stack"
  scale_fill_manual(values = c("Other"= "burlywood4", "Substrate"="bisque3", "Seagrass"="darkgreen", "Algae"="darkseagreen2", 
                               "CCA"="deeppink4", "Other_Inv"="midnightblue", "Sponge"="darksalmon", "Soft_Coral"="cyan3", 
                               "Hard_Coral"="deepskyblue4")) + 
  facet_grid(Site ~ ., scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change x-axis labels to percentage
  theme_bw() +
  labs(x = "Percentage", y = "Depth (m)") +
  scale_y_discrete(limits = rev(levels(df_long$Depth))) + # Depth from shallow to deep
  theme(axis.title.x = element_text(colour="Black", size=16),
        axis.title.y = element_text(colour="Black", size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=18), 
        legend.text = element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 18)
  )
# Print the plot
print(Groups_plot)


################ Creating DF of functional proportions groups per picture - 'Groups_all':

# to summarize how many live  and how many Non_live we have per picture:
Groups_all=aggregate(Count~SITE+DEPTH+PIC+GROUP, data=data_all, FUN="sum")

# to calculate the proportion:
Groups_all$prop=Groups_all$Count/64

# add variable that gives % coverage instead of proportion
Groups_all$Cov_percent=Groups_all$prop*100

# add variable that combines site and depth
Groups_all$site_depth=paste(Groups_all$SITE, Groups_all$DEPTH, sep="_")
str(Groups_all)

# Converting 'Groups_all' wide, adding 0 values and change it to long again:
Groups_all_image_wide=as.data.frame(Groups_all[,c(3,4,7)] %>%
                                       pivot_wider(names_from = GROUP, values_from = Cov_percent, values_fill=0))

Groups_all_image_wide$SITE=Groups_all$SITE[match(Groups_all_image_wide$PIC,Groups_all$PIC)]
Groups_all_image_wide$DEPTH=Groups_all$DEPTH[match(Groups_all_image_wide$PIC,Groups_all$PIC)]
str(Groups_all_image_wide)

Groups_all_image_wide_long=reshape2::melt(Groups_all_image_wide, id.vars=list("PIC", "SITE", "DEPTH"))

# Subset the 'Groups_all_image_wide_long' DF for the different live functional groups (7 different DF's):
str(Groups_all_image_wide_long)
HC_df <- subset(Groups_all_image_wide_long, variable == "Hard Coral")
SC_df <- subset(Groups_all_image_wide_long, variable == "Soft Coral")
Sponge_df <- subset(Groups_all_image_wide_long, variable == "Sponge")
Other_Inv_df <- subset(Groups_all_image_wide_long, variable == "Other_Inv")
Algae_df <- subset(Groups_all_image_wide_long, variable == "Algae")
CCA_df <- subset(Groups_all_image_wide_long, variable == "CCA")
Seagrass_df <- subset(Groups_all_image_wide_long, variable == "Seagrass")


################ Functional groups statistics & plotting (for each live group separately):
# 1. statistics: ART ANOVA test (Required packages: as above, lines: 172-173)
# 2. plotting box-plot % per depth (ggplot2)


####################### HARD CORALS: ART ANOVA test & plot:
### Aligned ranks anova:
model_HC = art(value ~ DEPTH + SITE + DEPTH:SITE, data = HC_df)
### Check the success of the procedure:
model_HC
### Conduct ANOVA:
anova(model_HC) # Significant for 3 factors!

# Posthoc pairwise comparison for SITE (if interaction significant, then posthoc for interaction is enough):
SITE_HC_posthoc = art.con(model_HC, "SITE")
SITE_HC_posthoc
Sum_SITE_HC = as.data.frame(SITE_HC_posthoc) # library(rcompanion)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_SITE_HC)

# Posthoc pairwise comparison for DEPTH (if interaction significant, then posthoc for interaction is enough)::
DEPTH_HC_posthoc = art.con(model_HC, "DEPTH")
DEPTH_HC_posthoc
Sum_DEPTH_HC = as.data.frame(DEPTH_HC_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_DEPTH_HC)

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_HC_posthoc = art.con(model_HC, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_HC_posthoc
### For Tukey-adjusted p-values, use adjust="tukey" otherwise use: "none"
Sum_DEPTH_SITE_HC = as.data.frame(DEPTH_SITE_HC_posthoc)
# Finding the statistical group letter: 
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_HC)


############### HC plot:
str(HC_df)

ggplot(HC_df, aes(x=DEPTH, y=value,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))+
  facet_grid(.~SITE) + 
  ggtitle("Hard Coral coverage (%)") +
  scale_x_discrete(limits = rev)+
  xlab("Depth (M)") + ylab("Precentage") +
  #scale_y_discrete(limits = rev(levels(HC_df$DEPTH))) + # Depth from shallow to deep
  theme_bw() +
  coord_flip()+
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


####################### SOFT CORALS: ART ANOVA test & plot:
model_SC = art(value ~ DEPTH + SITE + DEPTH:SITE, data = SC_df)
model_SC
anova(model_SC) # Significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_SC_posthoc = art.con(model_SC, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_SC_posthoc
Sum_DEPTH_SITE_SC = as.data.frame(DEPTH_SITE_SC_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_SC)

############### SC plot:

ggplot(SC_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("Soft Coral coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(SC_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


#### Soft corals LINE PLOT, library(ggplot2) & library(dplyr):
# Summarize the total counts for each depth
SC_totals_by_depth <- SC_labels_totals %>%
  group_by(DEPTH) %>%
  summarize(Total_Count = sum(Count))

# Create the point plot with a connecting line
ggplot(SC_totals_by_depth, aes(x = Total_Count, y = DEPTH)) +
  geom_point(size = 5, color = "blue") +  # Add points
  # geom_line(color = "blue") +  # Adding connecting line didn't work! Added in PP.
  scale_y_discrete(limits = rev(levels(SC_labels_totals$DEPTH))) +  # Reverse the y-axis for depth
  labs(x = "Total Number of Counts", y = "Depth (m)") +  # Labels
  theme_bw() +  # Use a minimal theme for clean visualization
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

#### Soft corals PIE CHART, library(ggplot2):
ggplot(SC_labels_totals, aes(x = "", y = Cov_percent, fill = LABEL_corrected)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "Lobophytum" = "burlywood4", "Nephthya" = "#b3de69", "Rhytisma" = "burlywood2", # d9d9d9 light grey
    "Sarcophyton" = "lemonchiffon2", "Sclerophytum" = "#80b1d3", "Xenia" = "deepskyblue4",
    "Rumphella" = "honeydew3", "SC" = "#bebada", "Tubipora" = "cyan3", 
    "Briareum" = "midnightblue", "Melithaea" = "#ccebc5", "Dendronephthya" = "#bc80bd")) + 
  theme_void() +  # Removes background, gridlines, and axis text
  labs(fill = "Soft coral") +  # Adds a legend title
  theme(legend.title = element_text(size=18), 
        legend.text = element_text(size = 14),  # Adjusts legend text size
  )


####################### CCA: ART ANOVA test & plot:
model_CCA = art(value ~ DEPTH + SITE + DEPTH:SITE, data = CCA_df)
model_CCA
anova(model_CCA) # Significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_CCA_posthoc = art.con(model_CCA, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_CCA_posthoc
Sum_DEPTH_SITE_CCA = as.data.frame(DEPTH_SITE_CCA_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_CCA)

############### CCA plot:

ggplot(CCA_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("CCA coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(CCA_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


####################### SPONGES: ART ANOVA test & plot:
model_Sponge = art(value ~ DEPTH + SITE + DEPTH:SITE, data = Sponge_df)
model_Sponge
anova(model_Sponge) # significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_Sponge_posthoc = art.con(model_Sponge, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_Sponge_posthoc
Sum_DEPTH_SITE_Sponge = as.data.frame(DEPTH_SITE_Sponge_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_Sponge)


############### Sponges plot:

ggplot(Sponge_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("Sponges coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(Sponge_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


####################### ALGAE: ART ANOVA test & plot:
model_Algae = art(value ~ DEPTH + SITE + DEPTH:SITE, data = Algae_df)
model_Algae
anova(model_Algae) # Significant for 2 factors: depth & interaction

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_Algae_posthoc = art.con(model_Algae, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_Algae_posthoc
Sum_DEPTH_SITE_Algae = as.data.frame(DEPTH_SITE_Algae_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_Algae)


############### Algae plot:

ggplot(Algae_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("Algae coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(Algae_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


####################### OTHER_INV: ART ANOVA test & plot:
model_Other_Inv = art(value ~ DEPTH + SITE + DEPTH:SITE, data = Other_Inv_df)
model_Other_Inv
anova(model_Other_Inv) # significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_Other_Inv_posthoc = art.con(model_Other_Inv, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_Other_Inv_posthoc
Sum_DEPTH_SITE_Other_Inv = as.data.frame(DEPTH_SITE_Other_Inv_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_Other_Inv)

############### Other_Inv plot:

ggplot(Other_Inv_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("Other_Inv coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(Other_Inv_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


####################### SEAGRASS ART ANOVA test & plot:
model_Seagrass = art(value ~ DEPTH + SITE + DEPTH:SITE, data = Seagrass_df)
model_Seagrass
anova(model_Seagrass) 

# Posthoc pairwise comparison for the interaction:
DEPTH_SITE_Seagrass_posthoc = art.con(model_Seagrass, "DEPTH:SITE", adjust="tukey")
DEPTH_SITE_Seagrass_posthoc
Sum_DEPTH_SITE_Seagrass = as.data.frame(DEPTH_SITE_Seagrass_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SITE_Seagrass)

############### Seagrass plot:

ggplot(Seagrass_df, aes(x=value, y=DEPTH,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + 
  ggtitle("Seagrass coverage (%)") +
  # scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Precentage") + ylab("Depth (M)") +
  scale_y_discrete(limits = rev(levels(Seagrass_df$DEPTH))) + # Depth from shallow to deep
  labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 


################# HARD CORALS 100% stacked bar plot:

###### Loading data (Tamar @ home desktop): 
setwd("C:/Users/Lenovo/Documents/STUDY/Data_Analayzing/Data_CSVs")
Hard_corals_plot<-read.csv('Hard_corals_plot.csv',header=TRUE, stringsAsFactors = TRUE)
str(Hard_corals_plot) 

######## Hard Corals plotting using library(ggplot2), library(reshape2):

# because we do not really want to tread depth as continuous value, 
# we will change it to be treated as factor:
Hard_corals_plot$Depth=as.factor(Hard_corals_plot$Depth)

# Reorder Site levels to "JG" at the top, "IUI-N" in the middle, "IUI-S" at the bottom
Hard_corals_plot$Site <- factor(Hard_corals_plot$Site, levels = c("JG", "IUI-N", "IUI-S"))

# Melt the data frame to long format for ggplot
Hard_corals_plot_long <- melt(Hard_corals_plot, id.vars = c("Site", "Depth"), variable.name = "Category", value.name = "Percentage")

# To change HC colors on plot using: library(RColorBrewer)
# Generate 44 distinct colors using the (RColorBrewer) package
# library(RColorBrewer)
num_colors <- 44
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)

# Reverse the order of the Category factor (the corals in the plot):
# library(forcats)
# HC_df_long$Category <- fct_rev(HC_df_long$Category)

# Create 100% stacked bar plot:
HC_plot <- ggplot(Hard_corals_plot_long, aes(x = Percentage, y = Depth, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) + 
  facet_grid(Site ~ ., scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +  # Change x-axis labels to percentage
  theme_bw() +
  labs(x = "Percentage", y = "Depth (m)", fill = "Hard Corals") + # Add axes and legend titles
  scale_y_discrete(limits = rev(levels(df_long$Depth))) + # Depth from shallow to deep
  scale_fill_manual(values = color_palette) +  # Use the generated color palette  
  theme(axis.title.x = element_text(colour="Black", size=16),
        axis.title.y = element_text(colour="Black", size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=16), 
        legend.text = element_text(size=10),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 18)) +
  guides(fill = guide_legend(ncol = 2))  # Arrange legend items in one or two columns
print(HC_plot)


############# HARD CORALS composition:

### prepare data frame with only counts of hard corals # NOTE: data_all has here already the column "LABEL_corrected"
str(data_all)

Hard_corals_all=data_all[data_all$GROUP=="Hard Coral",-c(6)]
str(Hard_corals_all)

Hard_corals_image=aggregate(Count~PIC+SITE+DEPTH+LABEL_corrected, data=Hard_corals_all, FUN="sum")
str(Hard_corals_image)

# To download to computer as csv: (row.names=FALSE is to cancel 1st column of numbers)
write.csv(Hard_corals_image, "Hard_corals_image.csv", row.names = FALSE)

Hard_corals_image_wide=as.data.frame(Hard_corals_image[,c(1,4,5)] %>%
                                       pivot_wider(names_from = LABEL_corrected, values_from = Count, values_fill=0))
  
row.names(Hard_corals_image_wide)=Hard_corals_image_wide$PIC
Hard_corals_image_wide=as.data.frame(t(Hard_corals_image_wide[,-1]))
str(Hard_corals_image_wide)

### count no. of observations of Hard corals per image and add as column: HC_counts to metadata DF:
# First, creating vector "ord" of all pics names: 
ord=as.vector(metadata$PIC)
Hard_corals_image_wide=Hard_corals_image_wide[,ord]

metadata$Site_Depth=paste(metadata$SITE, metadata$DEPTH, sep="_")
metadata$HC_counts=colSums(Hard_corals_image_wide)
str(metadata)

### Density plot of HC by sites:
ggplot(data = metadata, aes(x = HC_counts, fill = Site_Depth)) +
  geom_density(alpha = 0.7) + # Adjust alpha for transparency if needed
  #scale_fill_manual(values = c("JG" = "slategray", "IUI_S"="lightseagreen", "IUI_N" = "navyblue")) + 
  theme_bw() +
  #ylim(0,3.8) + # change limit of y axis
  # ggtitle("Density plot of live proportion by site") +
  xlab("HC_observations_per_image") + ylab("Density") +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        legend.position = c(1,1), # legend position in the top right of plot area
        legend.justification = c(1,1), # legend position inside the plot area
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  facet_grid(.~SITE)


########### Box plot of HC by site_depth:
ggplot(metadata, aes(x=Site_Depth, y=HC_counts))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))


########### Hard corals beta diversity statistics:

# Hellinger normalization:
Hard_corals_image_wide_hell=decostand(Hard_corals_image_wide, method="hellinger", MARGIN=2)

# Test contribution to variance by PERMANOVA
adonis2(t(Hard_corals_image_wide_hell)~DEPTH+SITE+DEPTH*SITE, data=metadata, distance="jaccard", permutations = 1000)

# Pairwise comparison:
pairwise.adonis2(t(Hard_corals_image_wide)~SITE, data=metadata, distance="jaccard")
pairwise.adonis2(t(Hard_corals_image_wide)~DEPTH, data=metadata, distance="jaccard")
pairwise.adonis2(t(Hard_corals_image_wide)~Site_Depth, data=metadata, distance="jaccard")


############# NMDS Calculate & plot for hard corals:

HC_NMDS=metaMDS(t(Hard_corals_image_wide_hell), k=3, try=100, trymax=100, distance="jaccard", autotransform=FALSE)
HC_NMDS$stress

metadata_HC_nmds=cbind(metadata, HC_NMDS$points)

ggplot(metadata_HC_nmds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=SITE, shape=SITE))

ggplot(metadata_HC_nmds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=DEPTH, shape=SITE))+
  scale_color_viridis(option="H", discrete=TRUE)+ # "H"-colors option for virdis, before it was "D"
  facet_grid(.~SITE)
# Need to run library(viridis) before.

# Plot HC_NMDS_1_2
ggplot(metadata_HC_nmds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=DEPTH, shape=SITE), size=5)+
  scale_color_viridis(option="H", discrete=TRUE)+
  scale_shape_manual(values=c(17,8,1)) +
theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
  ) 

# Plot HC_NMDS_1_3
ggplot(metadata_HC_nmds, aes(x=MDS1, y=MDS3))+
  geom_point(aes(color=DEPTH, shape=SITE), size=5)+
  scale_color_viridis(option="H", discrete=TRUE)+
  scale_shape_manual(values=c(17,8,1)) +
theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
  ) 


##### calculate LDM 
# install.packages("LDM")
library(LDM)

str(metadata)
library(dplyr)
library(stringr)
metadata <-  metadata %>% 
  mutate(start_location = str_replace(start_location, "-", "_"))


Hard_corals_image_wide_hell_t=as.data.frame(t(Hard_corals_image_wide_hell))

colSums(Hard_corals_image_wide_hell_t > 0)

Hard_corals_image_wide_hell_t_filter=Hard_corals_image_wide_hell_t[,colSums(Hard_corals_image_wide_hell_t > 0)>10]


HC_LDM=ldm(formula=Hard_corals_image_wide_hell_t_filter~DEPTH+SITE, 
           data=metadata, seed=12345, test.omni3 = TRUE, verbose=TRUE)

Hard_corals_covariates_global=as.data.frame(cbind(HC_LDM$p.global.freq,HC_LDM$p.global.tran, HC_LDM$p.global.pa,
                                                  HC_LDM$p.global.omni, HC_LDM$p.global.omni3))
colnames(Hard_corals_covariates_global)=c("freq","tran","pa","omni","omni3")

row.names(Hard_corals_covariates_global)=c("DEPTH","SITE")

HC_LDM_res= as.data.frame(rbind(HC_LDM$p.otu.omni3, HC_LDM$q.otu.omni3,
                                          HC_LDM$p.otu.omni, HC_LDM$q.otu.omni,
                                          HC_LDM$F.otu.pa,HC_LDM$p.otu.pa, HC_LDM$q.otu.pa,
                                          HC_LDM$F.otu.tran,HC_LDM$p.otu.tran, HC_LDM$q.otu.tran,
                                          HC_LDM$F.otu.freq, HC_LDM$p.otu.freq, HC_LDM$q.otu.freq))

vector_names=c("DEPTH_p.otu.omni3", "SITE_p.otu.omni3", 
               "DEPTH_q.otu.omni3","SITE_q.otu.omni3",
               "DEPTH_p.otu.omni","SITE_p.otu.omni",
               "DEPTH_q.otu.omni","SITE_q.otu.omni",
               "DEPTH_F.otu.pa","SITE_F.otu.pa",
               "DEPTH_p.otu.pa","SITE_p.otu.pa",
               "DEPTH_q.otu.pa","SITE_q.otu.pa",
               "DEPTH_F.otu.tran","SITE_F.otu.tran",
               "DEPTH_p.otu.tran","SITE_p.otu.tran",
               "DEPTH_q.otu.tran","SITE_q.otu.tran",
               "DEPTH_F.otu.freq","SITE_F.otu.freq",
               "DEPTH_p.otu.freq","SITE_p.otu.freq",
               "DEPTH_q.otu.freq","SITE_q.otu.freq")

row.names(HC_LDM_res)=vector_names

#MSRS_noldm_res$comparison=paste("MSRS_LDM")
#colnames(MSRS_noldm_res)=c("p.otu.omni3", "q.otu.omni3","p.otu.omni", "q.otu.omni","p.otu.pa", "q.otu.pa","p.otu.tran",
# "q.otu.tran","variable")
str(HC_LDM_res)
HC_LDM_res_t=as.data.frame(t(HC_LDM_res))

write.csv(HC_LDM_res_t,"LDM_HC_res.csv")


############### Code from Maya - plot HC after LDM test:

str(Hard_corals_image_wide)
Hard_corals_image_wide_1=Hard_corals_image_wide
Hard_corals_image_wide_1$coral=row.names(Hard_corals_image_wide_1)
Hard_corals_image_wide_1_melt=reshape2::melt(Hard_corals_image_wide_1, id.vars=list("coral"))
colnames(Hard_corals_image_wide_1_melt)=c("Coral", "PIC", "Prevalence")
Hard_corals_image_wide_1_melt$SITE=metadata$SITE[match(Hard_corals_image_wide_1_melt$PIC, metadata$PIC)]
Hard_corals_image_wide_1_melt$DEPTH=metadata$DEPTH[match(Hard_corals_image_wide_1_melt$PIC, metadata$PIC)]
Hard_corals_image_wide_1_melt$Site_Depth=metadata$Site_Depth[match(Hard_corals_image_wide_1_melt$PIC, metadata$PIC)]
str(Hard_corals_image_wide_1_melt)

write.csv(Hard_corals_image_wide_1_melt, "Hard_corals_image_wide_1_melt.csv")

keep_graphs=c("Acropora", "Montipora", "Porites", "Leptoseris", "Pocillopora", "Alveopora", "Stylophora")

str(Hard_corals_image_wide_1_melt)

Hard_corals_image_wide_1_melt_keep=Hard_corals_image_wide_1_melt[Hard_corals_image_wide_1_melt$Coral %in% keep_graphs,]
str(Hard_corals_image_wide_1_melt_keep)

ggplot(Hard_corals_image_wide_1_melt_keep, aes(x=Site_Depth, y=Prevalence))+
  geom_boxplot()+
  facet_wrap(Coral~., nrow=3, scales="free")
# Note: Wasn't used in the thesis, no need to run the code!


################ 7 Hard corals heat map: 

# Loading data of 7 selected HC for the heat map:
HC_Heatmap<-read.csv("HC_Heatmap.csv", header=TRUE, stringsAsFactors = TRUE)
setwd("C:/Users/Lenovo/Documents/STUDY/Data_Analayzing/R_scripts_etc")
HC_Heatmap$DEPTH=as.factor(HC_Heatmap$DEPTH) # setting depth as factor
HC_Heatmap$Site_Depth=paste(HC_Heatmap$SITE, HC_Heatmap$DEPTH, sep="_") # adding Site_Depth column
HC_Heatmap$SITE <- factor(HC_Heatmap$SITE, levels = c("IUI-S", "IUI-N", "JG")) # Reorder SITE factor levels

str(HC_Heatmap)

# Heat map of 7 corals with preferred order of corals & sites:
HC_heatmap <- ggplot(HC_Heatmap, aes(x=DEPTH, y=reorder(Coral,Order), fill=Prevalence)) +
  geom_tile() + 
  facet_grid(.~SITE, scales = "free", space = "free") +
scale_fill_viridis(option = "G") + # Use a viridis palette "G" for blue chosen colors
xlab("Depth (M)") + ylab("Corals") +
theme(axis.title.x = element_text(colour="Black", size=16),
      axis.title.y = element_text(colour="Black", size=16),
      axis.text.x = element_text(size=14),
      axis.text.y = element_text(size=14),
      legend.title = element_text(size=14), 
      legend.text = element_text(size=14),
      strip.text = element_text(size = 18) # Increase font size of facet labels
) 
print(HC_heatmap)

