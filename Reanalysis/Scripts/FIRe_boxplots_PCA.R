#### FIRE analysis
####### Eilat photo-physiology


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
library(multcompView)


# Loading data : 
setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")
FIRe_ave<- read.csv("FIRe_averaged_122025.csv", header=TRUE, stringsAsFactors=TRUE)
FIRe_ave$DEPTH <- as.factor(FIRe_ave$DEPTH)
str(FIRe_ave)

#### Removing all rows with NA's
FIRe_ave <- FIRe_ave %>% na.omit() 
str(FIRe_ave)


################ Normality tests for the fv_fm data: 


library(dplyr)
library(tidyr)
library(purrr)
library(broom)

str(FIRe_ave)

#### test normality using the Kolmogorov-Smirnov test
ks.normality <- FIRe_ave %>%
  pivot_longer(
    cols = 10:13,
    names_to = "variable",
    values_to = "value"
  ) %>%
  group_by(DEPTH,Species, variable) %>%
  summarise(
    ks = list(
      ks.test(
        value,
        "pnorm",
        mean(value, na.rm = TRUE),
        sd(value, na.rm = TRUE)
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(tidy_ks = purrr::map(ks, broom::tidy)) %>%
  unnest(tidy_ks) %>%
  dplyr::select(DEPTH,Species, variable,statistic, p.value) %>%
  dplyr::rename(
    D = statistic,
    P = p.value
  )

ks.normality= as.data.frame(ks.normality)


#### test homogeneity of variance

levene_results <- FIRe_ave %>%
  pivot_longer(cols = 10:13, names_to = "variable", values_to = "value") %>%
  group_by(Species, variable) %>%
  summarise(
    lv = list(car::leveneTest(value ~ as.factor(DEPTH))),
    .groups = "drop"
  ) %>%
  mutate(
    df1 = purrr::map_dbl(lv, ~ .x$Df[1]),
    df2 = purrr::map_dbl(lv, ~ .x$Df[2]),
    F   = purrr::map_dbl(lv, ~ .x$`F value`[1]),
    P   = purrr::map_dbl(lv, ~ .x$`Pr(>F)`[1])
  ) %>%
  dplyr::select(Species, variable, df1, df2, F, P)

levene_results
levene_results=as.data.frame(levene_results)


#### artanova tests 
# since for some we did not get normality and for others the variance was not homogeneous, we choose a non-parametric test for all.


FIRe_ave_long <- as.data.frame(FIRe_ave %>%
  pivot_longer(cols = 10:13, names_to = "variable", values_to = "responses"))
FIRe_ave_long$variable=factor(FIRe_ave_long$variable)
str(FIRe_ave_long)





art_results <- FIRe_ave_long %>%
  dplyr::mutate(
    Species = factor(Species),
    DEPTH   = factor(DEPTH),
    variable= factor(variable),
    Site    = factor(SITE)
  ) %>%
  group_by(Species, variable) %>%
  group_modify(~{
    dat <- drop_na(.x, responses, DEPTH, SITE)
    
    fit <- tryCatch(
      ARTool::art(responses ~ DEPTH + (1 | SITE), data = dat),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      return(tibble(
        term = NA_character_,
        df = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        note = fit$message
      ))
    }
    
    anova(fit) %>%
      as.data.frame() %>%
      rownames_to_column("term") %>%
      as_tibble() %>%
      dplyr::select(
        term,
        df        = any_of(c("Df", "df")),
        statistic = any_of(c("F", "F value", "F.value")),
        p.value   = any_of(c("Pr(>F)", "Pr(>Chisq)", "p.value"))
      ) %>%
      mutate(note = NA_character_)
  }) %>%
  ungroup()

art_results=as.data.frame(art_results)
art_results <- art_results[, !names(art_results) %in% c("note")]
write.csv(art_results, "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/FIRe_art_anova_results.csv")



######## posthoc art anova tests
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ARTool)

response_col <- "responses"   # <-- the numeric column in FIRe_ave_long
alpha <- 0.05
adj   <- "tukey"

make_form <- function(resp) {
  as.formula(paste0("`", resp, "` ~ DEPTH + (1|SITE)"))
}

fits <- FIRe_ave_long %>%
  mutate(
    Species  = factor(Species),
    variable = factor(variable),
    DEPTH    = factor(DEPTH),
    SITE     = factor(SITE)
  ) %>%
  group_by(Species, variable) %>%
  summarise(
    data = list(tidyr::drop_na(pick(everything()), all_of(c(response_col, "DEPTH", "SITE")))),
    .groups = "drop"
  ) %>%
  mutate(
    form = list(make_form(response_col)),
    art_fit = map(data, ~ tryCatch(ARTool::art(form[[1]], data = .x), error = identity))
  )

# omnibus p for DEPTH
omnibus <- fits %>%
  mutate(
    an = purrr::map(
      art_fit,
      function(x) {
        if (inherits(x, "error")) return(NULL)
        as.data.frame(anova(x))
      }
    ),
    p_depth = purrr::map_dbl(
      an,
      function(a) {
        if (is.null(a)) return(NA_real_)
        i <- which(rownames(a) == "DEPTH")
        if (length(i) != 1) return(NA_real_)
        pcol <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), colnames(a))[1]
        as.numeric(a[i, pcol])
      }
    ),
    note = purrr::map_chr(
      art_fit,
      function(x) {
        if (inherits(x, "error")) return(x$message)
        NA_character_
      }
    )
  ) %>%
  dplyr::select(Species, variable, p_depth, note)

posthoc <- fits %>%
  left_join(omnibus, by = c("Species", "variable")) %>%
  filter(!is.na(p_depth) & p_depth < alpha) %>%
  mutate(
    con = purrr::map(
      art_fit,
      function(x) ARTool::art.con(x, "DEPTH", adjust = adj)
    ),
    con_tbl = purrr::map(
      con,
      function(cc) tibble::as_tibble(summary(cc))
    )
  ) %>%
  dplyr::select(Species, variable, p_depth, con_tbl) %>%
  tidyr::unnest(con_tbl)


posthoc

write.csv(as.data.frame(posthoc), "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output/art_anova_posthoc_results.csv")


#posthoc contains the pairwise DEPTH contrasts only for the groups where dep_p < alpha.

##### now table of letters for the graph

p_col <- "p.value"   # change if your p-value column is named differently

cld_from_posthoc <- posthoc %>%
  dplyr::select(Species, variable, contrast, p = all_of(p_col)) %>%
  mutate(
    # normalize "A - B" / "A-B" -> "A-B"
    contrast = gsub("\\s+", "", as.character(contrast)),  # remove all spaces
    contrast = gsub("-", "-", contrast)                  # keep hyphen
  ) %>%
  group_by(Species, variable) %>%
  summarise(
    letters = list({
      df <- cur_data()
      
      # named p-value vector: names like "5-10"
      pvec <- df$p
      names(pvec) <- df$contrast
      
      # multcompLetters expects names "A-B" (no spaces)
      multcompView::multcompLetters(pvec)$Letters
    }),
    .groups = "drop"
  ) %>%
  unnest_longer(letters, indices_to = "DEPTH", values_to = "cld") %>%
  mutate(DEPTH = as.character(DEPTH)) %>%
  arrange(Species, variable, DEPTH)

cld_from_posthoc

cld_from_posthoc_wide=pivot_wider(cld_from_posthoc, names_from = DEPTH, values_from =cld )
write.csv(cld_from_posthoc_wide, "FIRe_art_test_posthoc_letters.csv")




##### boxplot graphs

str(FIRe_ave_long)
levels(FIRe_ave_long$Species)
FIRe_ave_long$Species=factor(FIRe_ave_long$Species, levels=c("Stylophora", "Porites", "Pocillopora", "Leptoseris"))

## remove outliers

FIRe_ave_long_no_out <- FIRe_ave_long %>%
  group_by(Species, variable) %>%   # drop "variable" if you don't have it
  mutate(
    q1 = quantile(responses, 0.25, na.rm = TRUE),
    q3 = quantile(responses, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lo = q1 - 2 * iqr,
    hi = q3 + 2 * iqr
  ) %>%
  filter(is.na(responses) | (responses >= lo & responses <= hi)) %>%  # keeps non-outliers
  ungroup() %>%
  dplyr::select(-q1, -q3, -iqr, -lo, -hi)


letters_df <- posthoc %>%
  dplyr::select(Species, variable, contrast, p = all_of(p_col)) %>%
  mutate(
    # normalize contrasts like "5 - 10" or "5-10" -> "5-10"
    contrast = gsub("\\s+", "", as.character(contrast))
  ) %>%
  group_by(Species, variable) %>%
  summarise(
    letters = list({
      pvec <- p
      names(pvec) <- contrast
      multcompView::multcompLetters(pvec)$Letters
    }),
    .groups = "drop"
  ) %>%
  unnest_longer(letters, indices_to = "DEPTH", values_to = "cld") %>%
  mutate(DEPTH = as.character(DEPTH))

normalize_depth <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x
}
str(letters_df)

letters_df$DEPTH_val <- gsub("DEPTH", "", letters_df$DEPTH)

label_df <- FIRe_ave_long_no_out %>%
  mutate(
    DEPTH_key = normalize_depth(DEPTH),
    Species   = as.character(Species),
    variable  = as.character(variable)
  ) %>%
  group_by(Species, variable, DEPTH_key) %>%
  summarise(y = max(responses, na.rm = TRUE), .groups = "drop") %>%
  group_by(Species, variable) %>%
  mutate(
    # offset per panel (works even with free_y)
    y = y + 0.05 * diff(range(y, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  left_join(
    letters_df %>%
      mutate(
        DEPTH_key = normalize_depth(as.numeric(DEPTH_val)),
        Species   = as.character(Species),
        variable  = as.character(variable)
      ) %>%
      dplyr::select(Species, variable, DEPTH_key, cld),
    by = c("Species", "variable", "DEPTH_key")
  )

##

FIRe_ave_long_no_out$Species=factor(FIRe_ave_long_no_out$Species, levels=c("Stylophora", "Porites", "Pocillopora", "Leptoseris"))
label_df$Species=factor(label_df$Species, levels=c("Stylophora", "Porites", "Pocillopora", "Leptoseris") )
letters_df$Species=factor(letters_df$Species, levels=c("Stylophora", "Porites", "Pocillopora", "Leptoseris") )

levels(FIRe_ave_long_no_out$variable)
FIRe_ave_long_no_out$variable=factor(FIRe_ave_long_no_out$variable, levels=c("fv_fm", "Pmax", "p", "Sigma"))
label_df$variable=factor(label_df$variable, levels=c("fv_fm", "Pmax", "p", "Sigma"))



library(ggh4x)

p_fire <- ggplot(FIRe_ave_long_no_out, aes(x = as.factor(DEPTH), y = responses)) +
  geom_boxplot(outlier.shape = NA, aes(fill=Species)) +
  scale_fill_manual(values=c("#aa87deff","#3d8abdff", "#a69c68ff", "#8aca6bff" ))+
  geom_text(
    data = label_df %>% filter(!is.na(cld)),
    aes(x = as.factor(DEPTH_key), y = y, label = cld),
    inherit.aes = FALSE,
    vjust = 0
  ) +
  facet_grid(variable~Species, scales = "free") +
  labs(x = "DEPTH", y = "responses")+
  facetted_pos_scales(
    y = list(
      variable == "fv_fm" ~ scale_y_continuous(limits = c(0, 0.8)),
      variable == "Pmax"  ~ scale_y_continuous(limits = c(0, 450)),
      variable == "p"    ~ scale_y_continuous(limits = c(0, 0.45)),
      variable == "Sigma" ~ scale_y_continuous(limits = c(100, 800))
    )
  )+
  theme_bw()+
  theme(axis.title.x = element_text(colour="Black", size=14),
        axis.title.y = element_text(colour="Black", size=14),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )

print(p_fire)
##
output_dir <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- Save as PNG ----
ggsave(
  filename = file.path(output_dir, "FIRe_photophysiology_by_depth.png"),
  plot = p_fire,
  width = 12,
  height = 8,
  dpi = 300
)

# ---- Save as PDF (journal-ready) ----
ggsave(
  filename = file.path(output_dir, "FIRe_photophysiology_by_depth.pdf"),
  plot = p_fire,
  width = 12,
  height = 8
)



#### fire_PCAs


library(FactoMineR)
library(factoextra)
library(viridis)
library(ggpubr)
library(gridExtra)

setwd("/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Data")
fire=read.csv("FIRe_averaged_122025_1.csv", header=TRUE,  stringsAsFactors = TRUE)
str(fire)
fire$ID=paste(fire$Species, fire$Colony, fire$DEPTH, fire$SITE, sep="_")
str(fire)
row.names(fire)=fire$ID

fire_1=fire[!is.na(fire$p),]


#Leptoseris
Leptoseris=fire_1[fire_1$Species=="Leptoseris",]
dim(Leptoseris)

Leptoseris.pca=PCA(Leptoseris[,c(10:13)], scale.unit = TRUE, ncp = 10, graph = TRUE)
Leptoseris_fviz_coord=as.data.frame(Leptoseris.pca$ind$coord)

Leptoseris_base_biplot <- fviz_pca_biplot(
  Leptoseris.pca,
  invisible = "ind", # This makes the individuals invisible
  label = "var",
  repel = TRUE,
  title = "Leptoseris",
  col.var = "black"
)

Leptoseris_final_biplot <- Leptoseris_base_biplot +
  geom_point(data=Leptoseris_fviz_coord, aes(x=Dim.1, y=Dim.2, size=as.factor(Leptoseris$of_live), alpha=as.factor(Leptoseris$of_live)), shape=19, color="#8aca6bff" ) +
  #scale_shape_manual(values=c(21,22,24))+
  scale_size_manual(values=c(1.5,2.5,3.5),labels = c("25m (0.7%)", "35m (2.6%)", "45m (2.9%)"))+
  scale_alpha_manual(values=c(1,0.7,0.5),labels = c("25m (0.7%)", "35m (2.6%)", "45m (2.9%)"))+
  scale_fill_viridis(option="G", discrete = TRUE, direction=-1)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme_minimal()+
  theme(legend.title = element_blank())

Leptoseris_final_biplot


#Pocillopora
Pocillopora=fire_1[fire_1$Species=="Pocillopora",]
dim(Pocillopora)

Pocillopora.pca=PCA(Pocillopora[,10:13], scale.unit = TRUE, ncp = 10, graph = TRUE)

Pocillopora_fviz_coord=as.data.frame(Pocillopora.pca$ind$coord)

Pocillopora_base_biplot <- fviz_pca_biplot(
  Pocillopora.pca,
  invisible = "ind", # This makes the individuals invisible
  label = "var",
  repel = TRUE,
  title = "Pocillopora",
  col.var="black"
)

Pocillopora_final_biplot <- Pocillopora_base_biplot +
  geom_point(data=Pocillopora_fviz_coord, aes(x=Dim.1, y=Dim.2, size=as.factor(Pocillopora$of_live), alpha=as.factor(Pocillopora$of_live)), shape=19, color="#a69c68ff" ) +
  #scale_shape_manual(values=c(21,22,24),)+
  scale_size_manual(values=c(1.5,2.5,3.5),labels = c("5m (1.0%)", "15m (2.3%)", "25m (2.9%)"))+
  scale_alpha_manual(values=c(1,0.7,0.5),labels = c("5m (1.0%)", "15m (2.3%)", "25m (2.9%)"))+
  scale_fill_viridis(option="G", discrete = TRUE, direction=-1)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme_minimal()+
  theme(legend.title = element_blank())

Pocillopora_final_biplot

# Porites
Porites=fire_1[fire_1$Species=="Porites",]
dim(Porites)

Porites.pca=PCA(Porites[,10:13], scale.unit = TRUE, ncp = 10, graph = TRUE)
print(Porites.pca)               
fviz_screeplot(Porites.pca, ncp=10)
plot.PCA(Porites.pca, axes = c(1,2), choix="ind")
fviz_pca_biplot(Porites.pca,  geom = "text")

Porites_fviz_coord=as.data.frame(Porites.pca$ind$coord)

Porites.pca$ind$coord

Porites_base_biplot <- fviz_pca_biplot(
  Porites.pca,
  invisible = "ind", # This makes the individuals invisible
  label = "var",
  repel = TRUE,
  title = "Porites",
  col.var="black"
)
str(Porites)
unique(log2(Porites$of_live))

Porites_final_biplot <- Porites_base_biplot +
  geom_point(data=Porites_fviz_coord, aes(x=Dim.1, y=Dim.2, size=as.factor(Porites$of_live), alpha=as.factor(Porites$of_live)), shape=19, color="#3d8abdff" ) +
  #scale_shape_manual(values=c(21,22,24))+
  scale_size_manual(values=c(1.5,2.2,2.5,3,3.5),labels = c("5m (0.34%)", "15m (2.5%)", "25m (3.6%)", "35m (5.1%)", "45m (7.7%)"))+
  scale_alpha_manual(values=c(1,0.8,0.6,0.5,0.4), labels = c("5m (0.34%)", "15m (2.5%)", "25m (3.6%)", "35m (5.1%)", "45m (7.7%)"))+
  scale_fill_viridis(option="G", discrete = TRUE, direction=-1)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme_minimal()+
  theme(legend.title = element_blank())

Porites_final_biplot

# Stylophora
Stylophora=fire_1[fire_1$Species=="Stylophora",]
dim(Stylophora)

Stylophora.pca=PCA(Stylophora[,10:13], scale.unit = TRUE, ncp = 10, graph = TRUE)
print(Stylophora.pca)               
fviz_screeplot(Stylophora.pca, ncp=10)
plot.PCA(Stylophora.pca, axes = c(1,2), choix="ind")
fviz_pca_biplot(Stylophora.pca,  geom = "text")

Stylophora_fviz_coord=as.data.frame(Stylophora.pca$ind$coord)

Stylophora_base_biplot <- fviz_pca_biplot(
  Stylophora.pca,
  invisible = "ind", # This makes the individuals invisible
  label = "var",
  repel = TRUE,
  title = "Stylophora",
  col.var="black"
)


str(Stylophora)
unique(Stylophora$of_live)

Stylophora_final_biplot <- Stylophora_base_biplot +
  geom_point(data=Stylophora_fviz_coord, aes(x=Dim.1, y=Dim.2, size=as.factor(Stylophora$of_live), alpha=as.factor(Stylophora$of_live)), shape=19, color="#aa87deff" ) +
  #scale_shape_manual(values=c(21,22,24))+
  scale_size_manual(values=c(1.5,2.2,2.5,3,3.5),labels = c("5m (1.6%)", "15m (2.3%)", "25m (2.9%)", "35m (6.8%)", "45m (7.4%)"))+
  scale_alpha_manual(values=c(1,0.8,0.6,0.5,0.4),labels = c("5m (1.6%)", "15m (2.3%)", "25m (2.9%)", "35m (6.8%)", "45m (7.4%)"))+
  scale_fill_viridis(option="G", discrete = TRUE, direction=-1)+
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  theme_minimal()+
  theme(legend.title = element_blank())

Stylophora_final_biplot

#######
library(cowplot)
library(ggplot2)

# ---- combine ----
combined_pca <- plot_grid(
  Stylophora_final_biplot, Porites_final_biplot, Pocillopora_final_biplot, Leptoseris_final_biplot,
  labels = c("A","B","C","D"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

# ---- output folder ----
output_dir <- "/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Coral-communities-along-the-depth-gradient-in-the-Red-Sea-reefs-of-Eilat/Reanalysis/Output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- reset devices ----
graphics.off()

# ---- save (try PNG) ----
ggsave(
  filename = file.path(output_dir, "PCA_combined_biplots.png"),
  plot = combined_pca,
  width = 10,
  height = 7,
  dpi = 300
)

# ---- save (PDF) ----
ggsave(
  filename = file.path(output_dir, "PCA_combined_biplots.pdf"),
  plot = combined_pca,
  width = 10,
  height = 7
)

# ---- open a new device and display ----
quartz()
print(combined_pca)
