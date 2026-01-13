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
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")




# Loading data 
setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")
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
p_livecov <- ggplot(live, aes(x=DEPTH, y=prop,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.1))+
  facet_grid(.~SITE) + # ggtitle("Live observations proportions in each site by depth") +
 scale_y_continuous(labels = percent_format(accuracy = 1)) +  # Change y-axis labels to percentage
  xlab("Depth (M)") + ylab("Precentage") +
  #labs(colour = "Depth (M)") + # Add the legend title
  theme_bw() +
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        #legend.title = element_text(size=14), 
        #legend.text = element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), # Remove minor grid lines
        strip.text = element_text(size = 18) # Increase font size of facet labels
  ) 
print(p_livecov)
# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Live_observations_by_depth_and_site.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Live_observations_by_depth_and_site.pdf"),
  width = 8,
  height = 6
)



############### Homogeneity of variance Levene's tests:

# install.packages("car")
library(car)
prop_data<-live$prop # on prop data of all sites & depths (just a trial) 

# For site groups:
#leveneTest(prop_data~SITE, data=live)

# For depth groups:
leveneTest(prop_data~DEPTH, data=live)


############## Density plots for live prop by site:

p_dens <-ggplot(data = live, aes(x = prop, fill = SITE)) +
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
print(p_dens)
# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Density plots for live prop by site.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Density plots for live prop by site.pdf"),
  width = 8,
  height = 6
)


###### live coral proportion differs among sites Kruskal–Wallis 



## ---- Libraries (quiet load optional) ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## You can either load FSA quietly OR not load it at all and just call FSA::dunnTest
suppressPackageStartupMessages(library(FSA))   # optional

## ---- Data checks ----
# Make sure SITE is a factor
live$SITE <- as.factor(live$SITE)

# Optional: drop missing values
live2 <- live %>%
  filter(!is.na(prop), !is.na(SITE))

# Sample sizes per site
print(table(live2$SITE))

## ---- 1) Overall test: Kruskal–Wallis ----
kw_site <- kruskal.test(prop ~ SITE, data = live2)
print(kw_site)

## ---- 2) Post hoc: Dunn test with multiple-comparison correction ----
# Choose method = "bh" (FDR) or "bonferroni"
dunn_site <- FSA::dunnTest(prop ~ SITE, data = live2, method = "bh")
print(dunn_site)

# Save Dunn results as a data frame (nice for reporting/export)
dunn_df <- dunn_site$res
print(dunn_df)

## ---- Optional: Save the stats table to CSV ----
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

write.csv(dunn_df,
          file = file.path(output_path, "Dunn_posthoc_SITE_live_prop.csv"),
          row.names = FALSE)



# Calculate median and quartiles
live_summary <- live %>%
  group_by(SITE) %>%
  summarise(
    n = n(),
    median_cover = median(prop, na.rm = TRUE),
    Q1 = quantile(prop, 0.25, na.rm = TRUE),
    Q3 = quantile(prop, 0.75, na.rm = TRUE)
  )

# View results
print(live_summary)

# Define output directory (same one you used before)
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Save to CSV
write.csv(
  live_summary,
  file = file.path(output_path, "LiveCoverage_median_quartiles_by_SITE.csv"),
  row.names = FALSE
)



########## Density plots for live prop by depth:

p_dens_depth <- ggplot(data = live, aes(x = prop, fill = DEPTH)) +
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
print(p_dens_depth)
# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Density plots for live prop by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Density plots for live prop by depth.pdf"),
  width = 8,
  height = 6
)
################ ART ANOVA test for the 'live' data:

# Required packages: install.packages("ARTool"), install.packages("emmeans"), install.packages("multcomp")
# and: install.packages("rcompanion"), install.packages("ggplot2"), install.packages("psych")

library(ARTool)

### Aligned ranks anova:
model_live <- art(prop ~ DEPTH, data = live)
anova(model_live)

### Check the success of the procedure:
model_live

### Conduct ANOVA:
anova(model_live)



# Posthoc pairwise comparison for DEPTH:
DEPTH_posthoc = art.con(model_live, "DEPTH")
DEPTH_posthoc
Sum_DEPTH = as.data.frame(DEPTH_posthoc)

# Finding the statistical group letter: 
live_depth_stat_letters=as.data.frame (cldList(p.value ~ contrast, data=Sum_DEPTH))

write.csv(Sum_DEPTH,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/live coverage depth posthoc.csv")
write.csv(live_depth_stat_letters,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/live coverage depth posthoc_letters.csv")




##### Composition data processing

setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")
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



######## Functional groups plotting using library(ggplot2), library(reshape2):
 
setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")
Fun_groups_plot<-read.csv('Fun_groups_plot.csv',header=TRUE, stringsAsFactors = TRUE)
str(Fun_groups_plot) 



# Change depth to factor:
Fun_groups_plot$Depth <- as.factor(Fun_groups_plot$Depth)

# Melt to long format (keep Site in the data, but we won't facet by it)
df_long <- melt(Fun_groups_plot,
                id.vars = c("Site", "Depth"),
                variable.name = "Category",
                value.name = "Percentage")

# Create 100% stacked bar plot (ONLY BY DEPTH)
p_groups_plot <- ggplot(df_long, aes(x = Percentage, y = Depth, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_fill_manual(values = c(
    "Other"      = "burlywood4",
    "Substrate"  = "bisque3",
    "Seagrass"   = "darkgreen",
    "Algae"      = "darkseagreen2",
    "CCA"        = "deeppink4",
    "Other_Inv"  = "midnightblue",
    "Sponge"     = "darksalmon",
    "Soft_Coral" = "cyan3",
    "Hard_Coral" = "deepskyblue4"
  )) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  labs(x = "Percentage", y = "Depth (m)") +
  scale_y_discrete(limits = rev(levels(df_long$Depth))) +
  theme(
    axis.title.x = element_text(colour = "Black", size = 16),
    axis.title.y = element_text(colour = "Black", size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 14)
  )

print(p_groups_plot)

# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Functional groups by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Functional groups by depth.pdf"),
  width = 8,
  height = 6
)

####### combing life coverage with function group

library(cowplot)

# Optional: remove legends if not needed
p_dens        <- p_dens        + theme(legend.position = "none")
p_dens_depth  <- p_dens_depth  + theme(legend.position = "none")
# keep legend for p_groups_plot if desired

# ---- Top row (A, B) ----
top_row <- plot_grid(
  p_dens,  # A
  p_dens_depth,        # B
  ncol = 2,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0, vjust = 1
)

# ---- Bottom row (C) ----
bottom_row <- plot_grid(
  p_groups_plot,
  labels = "C",
  label_size = 18,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0, vjust = 1
)

# ---- Combine into final panel ----
combined_panel <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 1.2)  # bottom panel slightly taller
)

# Preview
print(combined_panel)

output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

ggsave(
  filename = file.path(output_path, "Fig.1_combined.png"),
  plot = combined_panel,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_path, "Fig.1_combined.pdf"),
  plot = combined_panel,
  width = 10,
  height = 12
)



################ Creating DF of functional proportions groups per picture - 'Groups_all':

# to summarize how many live  and how many Non_live we have per picture:
Groups_all=aggregate(Count~DEPTH+PIC+GROUP, data=data_all, FUN="sum")

# to calculate the proportion:
Groups_all$prop=Groups_all$Count/64

# add variable that gives % coverage instead of proportion
Groups_all$Cov_percent=Groups_all$prop*100

# add variable that combines site and depth
Groups_all$site_depth=paste(Groups_all$SITE, Groups_all$DEPTH, sep="_")
str(Groups_all)

# Converting 'Groups_all' wide, adding 0 values and change it to long again:
Groups_all_image_wide=as.data.frame(Groups_all[,c(2,3,6)] %>%
                                       pivot_wider(names_from = GROUP, values_from = Cov_percent, values_fill=0))

#Groups_all_image_wide$SITE=Groups_all$SITE[match(Groups_all_image_wide$PIC,Groups_all$PIC)]
Groups_all_image_wide$DEPTH=Groups_all$DEPTH[match(Groups_all_image_wide$PIC,Groups_all$PIC)]
str(Groups_all_image_wide)

Groups_all_image_wide_long=reshape2::melt(Groups_all_image_wide, id.vars=list("PIC", "DEPTH"))

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
model_HC = art(value ~ DEPTH,  data = HC_df)
### Check the success of the procedure:
model_HC
### Conduct ANOVA:
anova(model_HC) # Significant for 3 factors!



# Posthoc pairwise comparison for DEPTH (if interaction significant, then posthoc for interaction is enough)::
DEPTH_HC_posthoc = art.con(model_HC, "DEPTH")
DEPTH_HC_posthoc
Sum_DEPTH_HC = as.data.frame(DEPTH_HC_posthoc)
# Finding the statistical group letter: 
HC_depth_stat_letters=cldList(p.value ~ contrast, data=Sum_DEPTH_HC)

write.csv(Sum_DEPTH_HC,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/HC_depth posthoc.csv")
write.csv(HC_depth_stat_letters,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/HC_depth posthoc_letters.csv")




############### HC plot:
str(HC_df)

p_hard <- ggplot(HC_df, aes(x=DEPTH, y=value,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))+
  #facet_grid(.~SITE) + 
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
print (p_hard)


# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Hard coral cover by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Hard coral cover by depth.pdf"),
  width = 8,
  height = 6
)

####################### SOFT CORALS: ART ANOVA test & plot:
model_SC = art(value ~ DEPTH,  data = SC_df)
model_SC
anova(model_SC) # Significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_SC_posthoc = art.con(model_SC, "DEPTH", adjust="tukey")
DEPTH_SC_posthoc
Sum_DEPTH_SC = as.data.frame(DEPTH_SC_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_SC)

SC_depth_stat_letters=cldList(p.value ~ contrast, data=Sum_DEPTH_SC)

write.csv(Sum_DEPTH_SC,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/SC_depth posthoc.csv")
write.csv(SC_depth_stat_letters,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/SC_depth posthoc_letters.csv")


############### SC plot:

p_soft <- ggplot(SC_df, aes(x=DEPTH, y=value,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))+
  #facet_grid(.~SITE) + 
  ggtitle("Soft Coral coverage (%)") +
  scale_x_discrete(limits = rev)+
  xlab("Depth (M)") + ylab("Precentage") +
  #scale_y_discrete(limits = rev(levels(SC_df$DEPTH))) + # Depth from shallow to deep
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
print(p_soft)



# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Soft coral cover by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Soft coral cover by depth.pdf"),
  width = 8,
  height = 6
)



####################### CCA: ART ANOVA test & plot:
model_CCA = art(value ~ DEPTH, data = CCA_df)
model_CCA
anova(model_CCA) # Significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_CCA_posthoc = art.con(model_CCA, "DEPTH", adjust="tukey")
DEPTH_CCA_posthoc
Sum_DEPTH_CCA = as.data.frame(DEPTH_CCA_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_CCA)

CCA_depth_stat_letters=cldList(p.value ~ contrast, data=Sum_DEPTH_CCA)

write.csv(Sum_DEPTH_CCA,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/CCA_depth posthoc.csv")
write.csv(CCA_depth_stat_letters,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/CCA_depth posthoc_letters.csv")



############### CCA plot:

p_cca <- ggplot(CCA_df, aes(x=DEPTH, y=value,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))+
  #facet_grid(.~SITE) + 
  ggtitle("CCA coverage (%)") +
  scale_x_discrete(limits = rev)+
  xlab("Depth (M)") + ylab("Precentage") +
  #scale_y_discrete(limits = rev(levels(CCA_df$DEPTH))) + # Depth from shallow to deep
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


print(p_cca)

# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "CCA cover by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "CCA cover by depth.pdf"),
  width = 8,
  height = 6
)
####################### SPONGES: ART ANOVA test & plot:
model_Sponge = art(value ~ DEPTH, data = Sponge_df)
model_Sponge
anova(model_Sponge) # significant for 3 factors!

# Posthoc pairwise comparison for the interaction:
DEPTH_Sponge_posthoc = art.con(model_Sponge, "DEPTH", adjust="tukey")
DEPTH_Sponge_posthoc
Sum_DEPTH_Sponge = as.data.frame(DEPTH_Sponge_posthoc)
cldList(p.value ~ contrast, data=Sum_DEPTH_Sponge)

Sponge_depth_stat_letters=cldList(p.value ~ contrast, data=Sum_DEPTH_Sponge)

write.csv(Sum_DEPTH_Sponge,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/Sponge_depth posthoc.csv")
write.csv(Sponge_depth_stat_letters,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/Sponge_depth posthoc_letters.csv")


############### Sponges plot:

p_sponge <-ggplot(Sponge_df, aes(x=DEPTH, y=value,colour=DEPTH))+
  geom_boxplot()+
  geom_point(position=position_jitter(width=0.2))+
  #facet_grid(.~SITE) + 
  ggtitle("Sponge coverage (%)") +
  scale_x_discrete(limits = rev)+
  xlab("Depth (M)") + ylab("Precentage") +
  #scale_y_discrete(limits = rev(levels(Sponge_df$DEPTH))) + # Depth from shallow to deep
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
print(p_sponge)
# Define output directory
output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
# Save plot 
ggsave(
  filename = file.path(output_path, "Spong cover by depth.png"),
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  filename = file.path(output_path, "Spong cover by depth.pdf"),
  width = 8,
  height = 6
)



########combine 4 plots to one

library(ggplot2)
library(cowplot)

# Remove legends (depth already on axis)

p_cca     <- p_cca     + theme(legend.position = "none")
p_hard    <- p_hard    + theme(legend.position = "none")
p_soft    <- p_soft    + theme(legend.position = "none")
p_sponge  <- p_sponge  + theme(legend.position = "none")

# Combine into one panel (2 columns x 3 rows)
combined_panel <- plot_grid(
    
  p_hard,   p_soft,
  p_sponge, p_cca,
  ncol = 2,
  align = "v", 
  labels = c("A", "B", "C", "D"),      # ⬅️ PANEL LABELS
  label_size = 18,                     # ⬅️ size of A–D
  label_fontface = "bold",             # ⬅️ bold labels
  label_x = 0.02,                      # ⬅️ left position
  label_y = 0.98,                      # ⬅️ top position
  hjust = 0, vjust = 1
)

# Preview
print(combined_panel)

output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

ggsave(
  filename = file.path(output_path, "BenthicCover_byDepth_combined.png"),
  plot = combined_panel,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_path, "BenthicCover_byDepth_combined.pdf"),
  plot = combined_panel,
  width = 10,
  height = 12
)


################# HARD CORALS 100% stacked bar plot:

##########HC by depth only

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(cowplot)

## ---- Data prep ----
Hard_corals_plot$Depth <- as.factor(Hard_corals_plot$Depth)

# Melt to long format
Hard_corals_plot_long <- melt(
  Hard_corals_plot,
  id.vars = c("Site", "Depth"),
  variable.name = "Category",
  value.name = "Percentage"
)

## ---- Aggregate across sites (mean per Depth x Category) ----
Hard_corals_depth <- aggregate(
  Percentage ~ Depth + Category,
  data = Hard_corals_plot_long,
  FUN = mean,
  na.rm = TRUE
)

## ---- Colors ----
num_colors <- length(unique(Hard_corals_depth$Category))
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)

## ---- Depth order (shallow to deep) ----
depth_levels <- levels(Hard_corals_plot_long$Depth)
Hard_corals_depth$Depth <- factor(Hard_corals_depth$Depth, levels = depth_levels)

## ---- Base plot (BY DEPTH ONLY) ----
HC_plot <- ggplot(Hard_corals_depth, aes(x = Percentage, y = Depth, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.6) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_reverse() +   # ⬅️ HIGH VALUES ON THE LEFT
  scale_y_discrete(limits = rev(depth_levels)) +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  labs(x = "Percentage", y = "Depth (m)", fill = NULL) +
  theme(
    axis.title.x = element_text(size = 12, colour = "black"),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  ) +
  guides(fill = guide_legend(ncol = 12, byrow = TRUE))

## ---- Extract legend robustly ----
leg_list <- get_plot_component(HC_plot, "guide-box", return_all = TRUE)

legend_grob <- NULL
for (g in leg_list) {
  if (!is.null(g) && inherits(g, "grob") &&
      !is.null(g$grobs) && length(g$grobs) > 0) {
    legend_grob <- g
    break
  }
}
if (is.null(legend_grob)) legend_grob <- leg_list[[length(leg_list)]]

legend_plot <- cowplot::ggdraw() + cowplot::draw_grob(legend_grob)

## ---- Main plot without legend ----
HC_plot_no_legend <- HC_plot + theme(legend.position = "none")

## ---- Optional: preview ----
print(HC_plot_no_legend)
grid.newpage(); grid.draw(legend_grob)   # preview legend grob
print(legend_plot)                       # preview legend plot


## ---- Output folder
output_dir <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


## ---- Save PDFs ----

# 1) Main plot: single journal column (narrow)

ggsave(
  filename = file.path(output_path, "HardCorals_main_single_column.png"),
  plot =HC_plot_no_legend,
  width = 5,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_path, "HardCorals_main_single_column.pdf"),
  plot = HC_plot_no_legend,
  width = 7,
  height = 5
)

# 2) Legend only: landscape (wide + short)


ggsave(
  filename = file.path(output_path, "HardCorals_legend_landscape.png"),
  plot =legend_grob,
  width = 12,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_path, "HardCorals_legend_landscape.pdf"),
  plot = legend_grob,
  width = 12,
  height = 5
)





#### Heat map of 7 corals with preferred order of corals & depth
library(ggplot2)
library(viridis)

# Load data
HC_Heatmap <- read.csv("HC_Heatmap.csv", header = TRUE, stringsAsFactors = TRUE)

# Set working directory (if needed)
setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")

# Set depth as factor
HC_Heatmap$DEPTH <- as.factor(HC_Heatmap$DEPTH)

# OPTIONAL: reorder DEPTH levels if needed (example: shallow to deep)
# HC_Heatmap$DEPTH <- factor(HC_Heatmap$DEPTH, levels = c("5","15","25","35","45"))

str(HC_Heatmap)

# ---- Heatmap by DEPTH only ----
HC_heatmap <- ggplot(
  HC_Heatmap,
  aes(x = DEPTH, y = reorder(Coral, Order), fill = Prevalence)
) +
  geom_tile() +
  scale_fill_viridis(option = "G") +
  xlab("Depth (m)") +
  ylab("Corals") +
  theme_bw() +
  theme(
    axis.title.x = element_text(colour = "Black", size = 16),
    axis.title.y = element_text(colour = "Black", size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14)
  )

print(HC_heatmap)

output_path <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

ggsave(
  filename = file.path(output_path, "HC_heatmap_depth.png"),
  plot = HC_heatmap,
  width = 10,
  height = 5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(output_path, "HC_heatmap_depth.pdf"),
  plot = HC_heatmap,
  width = 10,
  height = 5
)

                   






############# Create metadata table for images using library(vegan): 
str(data_all)
metadata=unique(data_all[,1:3])
#write.csv(metadata, "metadata.csv")

# Order the names (photos) in 2 DF's the same by creating the 'ord' vector:
row.names(metadata)=metadata$PIC
ord=match(colnames(data_all_image_wide), row.names(metadata))
metadata=metadata[ord,]



#### get_hart corals only
str(data_all)

Hard_corals_all=data_all[data_all$GROUP=="Hard Coral",-c(6)]
str(Hard_corals_all)

Hard_corals_image=aggregate(Count~PIC+SITE+DEPTH+LABEL_corrected, data=Hard_corals_all, FUN="sum")
str(Hard_corals_image)

Hard_corals_image_wide=as.data.frame(Hard_corals_image[,c(1,4,5)] %>%
                                       pivot_wider(names_from = LABEL_corrected, values_from = Count, values_fill=0))


row.names(Hard_corals_image_wide)=Hard_corals_image_wide$PIC
Hard_corals_image_wide=as.data.frame(t(Hard_corals_image_wide[,-1]))
str(Hard_corals_image_wide)

### count no. of observations of Hard corals per image and add as column to metadata

metadata$Site_Depth=paste(metadata$SITE, metadata$DEPTH, sep="_")
metadata$HC_counts=colSums(Hard_corals_image_wide)
str(metadata)



##### calculate LDM 
Hard_corals_image_wide_hell=decostand(Hard_corals_image_wide, method="hellinger", MARGIN=2)

library(LDM)

str(metadata)
library(dplyr)
library(stringr)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("BiocParallel")
#install.packages("LDM")   # or BiocManager::install("LDM") if needed
#library(LDM)


Hard_corals_image_wide_hell_t=as.data.frame(t(Hard_corals_image_wide_hell))

Hard_corals_image_wide_hell_t_filter=Hard_corals_image_wide_hell_t[,colSums(Hard_corals_image_wide_hell_t > 0)>10]


HC_LDM=ldm(formula=Hard_corals_image_wide_hell_t_filter~DEPTH, 
           data=metadata, seed=12345, test.omni3 = TRUE, verbose=TRUE)


Hard_corals_covariates_global=as.data.frame(cbind(HC_LDM$p.global.freq,HC_LDM$p.global.pa))

colnames(Hard_corals_covariates_global)=c("freq","pa")

row.names(Hard_corals_covariates_global)=c("DEPTH")


HC_LDM_res= as.data.frame(rbind(HC_LDM$F.otu.pa,HC_LDM$p.otu.pa, HC_LDM$q.otu.pa,
                                HC_LDM$F.otu.freq, HC_LDM$p.otu.freq, HC_LDM$q.otu.freq))


vector_names=c("DEPTH_F.otu.pa",
               "DEPTH_p.otu.pa",
               "DEPTH_q.otu.pa",
               "DEPTH_F.otu.freq",
               "DEPTH_p.otu.freq",
               "DEPTH_q.otu.freq")

row.names(HC_LDM_res)=vector_names

#MSRS_noldm_res$comparison=paste("MSRS_LDM")
#colnames(MSRS_noldm_res)=c("p.otu.omni3", "q.otu.omni3","p.otu.omni", "q.otu.omni","p.otu.pa", "q.otu.pa","p.otu.tran", "q.otu.tran","variable")
str(HC_LDM_res)
HC_LDM_res_t=as.data.frame(t(HC_LDM_res))


write.csv(HC_LDM_res_t,"LDM_HC_res.csv")


#### HC_relative abundances by depth (from total live coverage)

str(Hard_corals_image)
row.names(Hard_corals_image_wide)

Hard_corals_image_stats=as.data.frame(Hard_corals_image %>%
                                        group_by(LABEL_corrected, DEPTH) %>%
                                        summarise(count=sum(Count),
                                                  prevalence=sum(Count != 0)))


# total counts live by depth
str(data_all)

live_by_depth= data_all[data_all$Live_Non_live=="Live",] %>%
  group_by(DEPTH) %>%
  summarise(Total_live=sum(Count))

str(Hard_corals_image_stats)
Hard_corals_image_stats$Live_total=live_by_depth$Total_live[match(Hard_corals_image_stats$DEPTH,live_by_depth$DEPTH )]

Hard_corals_image_stats$Percent_of_Live=Hard_corals_image_stats$count/Hard_corals_image_stats$Live_total*100

Hard_corals_depth_percent_live=pivot_wider(Hard_corals_image_stats[,c(1,2,6)],
                                           names_from = "DEPTH",
                                           values_from = "Percent_of_Live",
                                           values_fill = list(value = 0))
str(Hard_corals_depth_percent_live)
colnames(Hard_corals_depth_percent_live)=c("LABEL_corrected", "m5", "m15", "m25", "m35", "m45")

HC_LDM_res_t$D5=Hard_corals_depth_percent_live$m5[match(row.names(HC_LDM_res_t),Hard_corals_depth_percent_live$LABEL_corrected )]
HC_LDM_res_t$D15=Hard_corals_depth_percent_live$m15[match(row.names(HC_LDM_res_t),Hard_corals_depth_percent_live$LABEL_corrected )]
HC_LDM_res_t$D25=Hard_corals_depth_percent_live$m25[match(row.names(HC_LDM_res_t),Hard_corals_depth_percent_live$LABEL_corrected )]
HC_LDM_res_t$D35=Hard_corals_depth_percent_live$m35[match(row.names(HC_LDM_res_t),Hard_corals_depth_percent_live$LABEL_corrected )]
HC_LDM_res_t$D45=Hard_corals_depth_percent_live$m45[match(row.names(HC_LDM_res_t),Hard_corals_depth_percent_live$LABEL_corrected )]

str(HC_LDM_res_t)


write.csv(HC_LDM_res_t,"/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/HC_LDM_res_t.csv")
