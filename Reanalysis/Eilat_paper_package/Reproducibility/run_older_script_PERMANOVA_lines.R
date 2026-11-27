suppressMessages({library(tidyverse);library(vegan)})
setwd("<folder with the data files>")
lc_sf.data<-read.csv('Annotations_MO_25062024_ver1.csv',header=TRUE, stringsAsFactors = TRUE)
str(lc_sf.data)
levels(lc_sf.data$GROUP)

# because we do not really want to tread depth as continuous value, 
# we will change it to be treated as factor:
lc_sf.data$DEPTH=as.factor(lc_sf.data$DEPTH)
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
Hard_corals_image_wide_hell=decostand(Hard_corals_image_wide, method="hellinger", MARGIN=2)
cat("CLASSES:\n"); print(sapply(metadata[,c("PIC","SITE","DEPTH","Site_Depth")],class)); print(table(metadata$DEPTH)); print(all(colnames(Hard_corals_image_wide)==rownames(metadata)))
X<-t(Hard_corals_image_wide_hell); mm<-model.matrix(~DEPTH+SITE+DEPTH*SITE,data=metadata); cat("model.matrix cols",ncol(mm),"rank",qr(mm)$rank,"\n")
print(adonis2(X~DEPTH+SITE+DEPTH*SITE, data=metadata, distance="jaccard", permutations = 1000))
print(adonis2(X~DEPTH+SITE+DEPTH*SITE, data=metadata, distance="jaccard", permutations = 1000, by="terms"))
a1<-adonis2(X~DEPTH*SITE,data=metadata,distance="jaccard",permutations=99);a2<-adonis2(X~DEPTH*SITE,data=metadata,method="bray",permutations=99);cat("SS distance arg",a1$SumOfSqs[1],"method bray",a2$SumOfSqs[1],"\n")
sessionInfo()
