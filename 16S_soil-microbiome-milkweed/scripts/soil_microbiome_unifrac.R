####Weighted Unifrac Distances for Soil Microbiome 6/25/2025
library(ape) 
library(vegan) 
library(ggplot2)
library(tidyverse)

#------Read in Unifrac Distance Matrix
weighted_unifrac <- read.csv("weighted-distance-matrix.csv", header = T, row.names = 1)

# First, make sure both have row names set correctly
# (optional: inspect the current row names)
head(rownames(weighted_unifrac))
head(rownames(metadata))

# Match metadata rows to genus_table rows (reorder metadata)
metadata_matched <- metadata[rownames(weighted_unifrac), ]

# Confirm the row names are now aligned
all(rownames(weighted_unifrac) == rownames(metadata_matched))

PCOA <- data.frame(pcoa(weighted_unifrac)$vectors)

new_names <- rep("", ncol(PCOA))

for(i in 1:ncol(PCOA)){ 
  new_names[i] <- paste("PC",i, sep="") 
}

names(PCOA) <- new_names

PCOA$SampleID <- rownames(PCOA)

metadata_matched$SampleID <- rownames(metadata_matched)

PCOA <- merge(PCOA, metadata_matched, by = "SampleID")


ggplot(PCOA) + geom_point(aes(x = PC1, y = PC2, color = Group)) + labs(title="PCoA Plot")

order <- factor(PCOA$Group, level = c("Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
                                      "Custer SY", "Custer SP", "Custer Control"))

my_beta_colors <- c(
  "Sioux Falls SY" = "darkblue",
  "Sioux Falls SP" = "turquoise",
  "Sioux Falls Control" = "#9ecae1",
  "Custer SY" = "#a50f15",
  "Custer SP" = "violet",
  "Custer Control" = "#fcae91"
)

png("Unifrac_PCoA_weighted.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(PCOA, aes(x=PC1, y=PC2, color=order, fill=order)) + geom_point(shape=21, color="black",size=4) + scale_fill_manual(values=my_beta_colors) + xlab(paste("PCo1 - ", pc1, "%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Unifrac Distances (weighted)") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(weighted_unifrac ~ metadata_matched$Group)
########### Df SumOfSqs      R2     F Pr(>F)    
##Model     5  0.10968 0.87457 46.02  0.001 ***
## Residual 33  0.01573 0.12543                 
###Total    38  0.12541 1.00000 


#------Read in Unifrac Distance Matrix
unweighted_unifrac <- read.csv("unweighted_unifrac_distance-matrix.csv", header = T, row.names = 1)

# First, make sure both have row names set correctly
# (optional: inspect the current row names)
head(rownames(unweighted_unifrac))
head(rownames(metadata))

# Match metadata rows to genus_table rows (reorder metadata)
metadata_matched <- metadata[rownames(unweighted_unifrac), ]

# Confirm the row names are now aligned
all(rownames(unweighted_unifrac) == rownames(metadata_matched))

PCOA2 <- data.frame(pcoa(unweighted_unifrac)$vectors)

new_names <- rep("", ncol(PCOA2))

for(i in 1:ncol(PCOA2)){ 
  new_names[i] <- paste("PC",i, sep="") 
}

names(PCOA2) <- new_names

PCOA2$SampleID <- rownames(PCOA2)

metadata_matched$SampleID <- rownames(metadata_matched)

PCOA2 <- merge(PCOA2, metadata_matched, by = "SampleID")


ggplot(PCOA2) + geom_point(aes(x = PC1, y = PC2, color = Group)) + labs(title="PCoA Plot")

order <- factor(PCOA2$Group, level = c("Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
                                       "Custer SY", "Custer SP", "Custer Control"))

png("Unifrac_PCoA_UNweighted.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(PCOA2, aes(x=PC1, y=PC2, color=order, fill=order)) + geom_point(shape=21, color="black",size=4) + scale_fill_manual(values=my_beta_colors) + xlab(paste("PCo1 - ", pc1, "%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Unifrac Distances (unweighted)") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(unweighted_unifrac ~ metadata_matched$Group )
############# Df SumOfSqs      R2      F Pr(>F)    
#####Model     5   2.5195 0.23763 2.0572  0.001 ***
  ##Residual 33   8.0832 0.76237                  
###Total    38  10.6027 1.00000 


##########################Adding Alpha Diversity - Faith PD for fun
library(ape) 
library(vegan) 
library(ggplot2)

#------Read in Unifrac Distance Matrix
faith <- read.csv("faith-alpha-diversity.csv", header = T, row.names = 1)

# First, make sure both have row names set correctly
# (optional: inspect the current row names)
head(rownames(faith))
head(rownames(metadata))

# Match metadata rows to genus_table rows (reorder metadata)
metadata_matched <- metadata[rownames(faith), ]

# Confirm the row names are now aligned
all(rownames(faith) == rownames(metadata_matched))

names(faith) <- new_names

faith$SampleID <- rownames(faith)

metadata_matched$SampleID <- rownames(metadata_matched)


#####merge data frames for ggplot
faith_pd <- merge(faith, metadata_matched, by = "SampleID")

###Plot
level_order <- factor(faith_pd$Group, level = c("Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
                                                "Custer SY", "Custer SP", "Custer Control"))

faith_pd %>% 
  filter(Group %in% c("Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
                      "Custer SY", "Custer SP", "Custer Control")) %>%
  ggplot(aes(x=level_order, y=faith_pd, fill=level_order)) +
  geom_boxplot() 

my_colors <- c(
  "Sioux Falls SY" = "darkblue",
  "Sioux Falls SP" = "#5dade2",
  "Sioux Falls Control" = "#aed6f1",
  "Custer SY" = "#e74c3c",
  "Custer SP" = "#f1948a",
  "Custer Control" = "#fadbd8")


png("faith_pd.png", height = 4, width = 6, units = 'in', res = 600)
faith_pd %>% 
  filter(Group %in% c("Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
                      "Custer SY", "Custer SP", "Custer Control")) %>%
  ggplot(aes(x=level_order, y=faith_pd, fill=factor(level_order))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) + ggtitle("Alpha Diversity (Faith's Phylogenetic Diversity)") + labs(x="Group",y="Faith's PD") + theme_classic() + scale_color_manual(values=my_colors) + scale_fill_manual(values=my_colors) + theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = 7, colour = "black"), axis.text.y = element_text(size = 7, colour = "black"))
dev.off()

shapiro.test(faith_pd$faith_pd)
###p-value = 0.008

kruskal.test(faith_pd$faith_pd~faith_pd$Group)
##p-value=0.05


faith_pd$Group <- factor(faith_pd$Group, levels = c(
  "Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
  "Custer SY", "Custer SP", "Custer Control"
))

png("Faith_genus_withstars.png", height = 5, width = 6, units = 'in', res = 600)
ggboxplot(faith_pd, x = "Group", y = "faith_pd", fill = "Group", palette = my_colors) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 4,
                     comparisons = list(
                       c("Sioux Falls SY", "Custer SY"),
                       c("Sioux Falls SP", "Custer SP"),
                       c("Custer SP", "Custer Control")
                     ),
                     label.y = c(262.40, 262.45, 262.50)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(title = "Faith's PD", y = "Faith PD", x = "Group")
dev.off()
