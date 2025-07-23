######Soil Microbiome - ITS -Genus Level Analysis - ITS 07/15/25 (with removed samples O6, P4,)

library (vegan)
library (ape)
library (pgirmess)
library (labdsv)
library (ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ColorBrewer)

#########################################################################################
#----  ITS Alpha-diversity Analysis 
#########################################################################################

#------Read in Genus Table table
genus_table_in <- read.csv("ITS_genus_table.csv", header = T, row.names = 1)

#-- Additional steps to filter
genus.sums <- colSums(genus_table_in)
genus_filtered <- genus_table_in[ , which(genus.sums > 5)] #---Check ASV sum should be more than 5
genus_filtered_filtered <- dropspc(genus_filtered, 2) #--- ASV should be present in more than 2 samples
genus_table <- genus_filtered_filtered[order(rownames(genus_filtered_filtered)),]

#---- read in metadata -
metadata <- read.csv("soil_metadata_clean.txt", sep="\t", row.names = 1)

###merge metadata columns to break down soil type and milkweed species
metadata$Group = paste(metadata$soil.type, metadata$species, sep=" ")

# First, make sure both have row names set correctly
# (optional: inspect the current row names)
head(rownames(genus_table))
head(rownames(metadata))

# Match metadata rows to genus_table rows (reorder metadata)
metadata_matched <- metadata[rownames(genus_table), ]

# Confirm the row names are now aligned
all(rownames(genus_table) == rownames(metadata_matched))

###alpha diversity
Obs <- rowSums(genus_table > 0)
shann <- diversity(genus_table)
InvSimp <- diversity(genus_table, "invsimpson")
Obs.rare<- rarefy(genus_table, min(rowSums(genus_table)))
total_diversity <- data.frame(Obs,Obs.rare,shann,InvSimp)

total_diversity_new <- merge(total_diversity, metadata_matched, by=0, all=F)
rownames(total_diversity_new) <- total_diversity_new$Row.names; total_diversity_new$Row.names <- NULL

####################################
#----ALPHA DIVERSITY
####################################
Observed <- total_diversity_new[,c(1,7)]

library(dplyr) 

#convert Group column to factor
level_order <- factor(Observed$Group, level = c('Sioux Falls SY', 'Sioux Falls SP', 'Sioux Falls Control', 'Custer SY', 'Custer SP', 'Custer Control'))

Observed %>% 
  filter(Group %in% c('Sioux Falls SY','Sioux Falls SP','Sioux Falls Control','Custer SY','Custer SP','Custer Control')) %>%
  ggplot(aes(x=level_order, y=Obs, fill=level_order)) +
  geom_boxplot() 

ITS_colors <- c("#1b9e77",  # forest green
                "#d95f02",  # earthy orange
                "#7570b3",  # deep purple
                "#e7298a",  # pink-fuchsia
                "#66a61e",  # olive green
                "#e6ab02")  # golden brown

png("ITS_Observed.png", height = 4, width = 6, units = 'in', res = 600)
Observed %>% 
  filter(Group %in% c('Sioux Falls SY','Custer SY','Sioux Falls SP','Custer SP','Custer Control','Sioux Falls Control')) %>%
  ggplot(aes(x=level_order, y=Obs, fill=level_order, varwidth=TRUE)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) + ggtitle("Observed genera") + labs(x="Group", y="Observed genera") + theme_classic() + scale_color_manual(values=ITS_colors) + scale_fill_manual(values=ITS_colors) + theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = 6, colour = "black"), axis.text.y = element_text(size = 7, colour = "black"))
dev.off()

shapiro.test(Observed$Obs)
##p-value=0.002

kruskal.test(Obs ~ Group, data = Observed)
#P-value=6.972e-06

pairwise.wilcox.test(Observed$Obs, Observed$Group,
                     p.adjust.method = "BH", exact = FALSE)  # Benjamini-Hochberg for multiple testing

###Try a different plot for fun
library(ggpubr)
library(ggplot2)

Observed$Group <- factor(Observed$Group, levels = c(
  "Sioux Falls SY",
  "Sioux Falls SP",
  "Sioux Falls Control",
  "Custer SY",
  "Custer SP",
  "Custer Control"
))


png("ITS_Observed_genera_pvalues_withstars.png", height = 5, width = 6, units = 'in', res = 600)
ggboxplot(Observed, x = "Group", y = "Obs", fill = "Group", palette = ITS_colors) +
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("Sioux Falls SP", "Custer SP"),
                       c("Sioux Falls SY", "Custer SY"),
                       c("Sioux Falls Control", "Custer Control")
                     ),
                     p.adjust.method = "BH",
                     label = "p.signif",
                     size = 4) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "right") +
  ggtitle("Observed Genera")
dev.off()

####Shannon
Shannon <- total_diversity_new[,c(3,7)]
level_order <- factor(Shannon$Group, level = c('Sioux Falls SY','Sioux Falls SP','Sioux Falls Control','Custer SY','Custer SP','Custer Control'))

Shannon$Group <- factor(Shannon$Group, levels = c(
  "Sioux Falls SY",
  "Sioux Falls SP",
  "Sioux Falls Control",
  "Custer SY",
  "Custer SP",
  "Custer Control"
))


###Shannon without log2
png("ITS_Shannon_genus.png", height = 4, width = 6, units = 'in', res = 600)
Shannon %>% 
  filter(Group %in% c('Sioux Falls SY','Sioux Falls SP','Sioux Falls Control','Custer SY','Custer SP','Custer Control')) %>%
  ggplot(aes(x=level_order, y=shann, varwidth=TRUE, fill=factor(level_order))) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) + ggtitle("Alpha Diversity (Shannon Index)") + labs(x="Group",y="Shannon Index") + theme_classic() + scale_color_manual(values=ITS_colors) + scale_fill_manual(values=ITS_colors) + theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = 7, colour = "black"), axis.text.y = element_text(size = 7, colour = "black"))
dev.off()

shapiro.test(Shannon$shann)
###p-value = 0.01

kruskal.test(shann ~ Group, data = Shannon)
#p-value =5.035e-05

pairwise.wilcox.test(Shannon$shann, Shannon$Group,
                     p.adjust.method = "BH", exact = FALSE)  # Benjamini-Hochberg for multiple testing


# Ensure Group is a factor in the right order
Shannon$Group <- factor(Shannon$Group, levels = c(
  "Sioux Falls SY", "Sioux Falls SP", "Sioux Falls Control",
  "Custer SY", "Custer SP", "Custer Control"
))



# Boxplot with ANOVA and Tukey HSD (stars)

png("ITS_Shannon_genus_withstars.png", height = 5, width = 6, units = 'in', res = 600)
ggboxplot(Shannon, x = "Group", y = "shann", fill = "Group", palette = ITS_colors) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(
                       c("Sioux Falls SY", "Custer SY"),
                       c("Sioux Falls SP", "Custer SP"),
                       c("Sioux Falls Control", "Custer Control")
                     ),
                     p.adjust.method = "BH",
                     label = "p.signif",
                     size = 4) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "right") +
  ggtitle("Shannon Diversity Index")
dev.off()

############################################
#-- Beta diversity PCoA plot
############################################
genus_table_2 <- merge(genus_table, metadata_matched, by=0, all=F)
rownames(genus_table_2) <- genus_table_2$Row.names; genus_table_2$Row.names <- NULL

genera <- genus_table_2[,1:114]
meta <- genus_table_2[,115:117]

genera_relab <- decostand(genera, method = "total")*100
write.table (genera_relab, file = "ITS_relative_adundance_GENUS.txt", sep = "\t")
bray_dist <- vegdist(genera_relab, method = "bray", binary = TRUE)
bray_pcoa <- pcoa (bray_dist)
bray_pcoa$values[1:2,]
pc1 <- round(bray_pcoa$values$Rel_corr_eig[1]*100, 2)
pc2 <- round(bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, meta)

##Shapes <- c(16, 17, 18, 15)

order <- factor(Bray_PCoA_MATRIX_New$Group, level = c('Sioux Falls SY', 'Sioux Falls SP', 'Sioux Falls Control', 'Custer SY', 'Custer SP', 'Custer Control'))

png("ITS_Bray_PCoA_Unweighted_genus.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, color=order, fill=order )) + geom_point(shape=21, color="black",size=4) + scale_fill_manual(values=ITS_colors) + xlab(paste("PCo1 - ", pc1,"%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances (unweighted)") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(formula = bray_dist ~ meta$Group)
########Df SumOfSqs      R2      F Pr(>F)    
#Model     5  0.74161 0.80587 28.228  0.001 ***
##Residual 34  0.17865 0.19413                  
##Total    39  0.92026 1.00000                  
#


####################Bray-Curtis Weighted
bray_dist <- vegdist(genera_relab, method = "bray", binary = FALSE)
bray_pcoa <- pcoa (bray_dist)
bray_pcoa$values[1:2,]
pc1 <- round(bray_pcoa$values$Rel_corr_eig[1]*100, 2)
pc2 <- round(bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, meta)


png("ITS_Bray_PCoA_WEIGHTED.png", height = 4, width = 6, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, color=order, fill=order)) + geom_point(shape=21, color="black",size=4) + scale_fill_manual(values=ITS_colors) + xlab(paste("PCo1 - ", pc1, "%", sep="")) + ylab(paste("PCo2 - ", pc2, "%", sep="")) + ggtitle("Bray-Curtis Distances (weighted)") + 
  theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

adonis2(bray_dist~meta$Group)
###adonis2(formula = bray_dist ~ meta$Group)
#######Df SumOfSqs     R2     F Pr(>F)    
#Model     5   1.3589 0.7725 23.09  0.001 ***
##Residual 34   0.4002 0.2275                 
##Total    39   1.7591 1.0000 

  # Observed Richness
richness <- specnumber(genus_table) 

### Shannon's H'

H <- diversity(genus_table)

# Pielou's Evenness
evenness <- H/log(richness)

# Create alpha diversity dataframe including environmental data
alpha <- cbind(shannon = H, richness = richness, pielou = evenness, meta)
head(alpha)

total_alpha_new <- merge(alpha, metadata_matched, by=0, all=F)
rownames(total_alpha_new) <- total_alpha_new$Row.names; total_alpha_new$Row.names <- NULL


####Pielou's evenness
Pielou <- total_alpha_new[,c(3,6)]
level_order <- factor(Pielou$Group.x, level = c('Sioux Falls SY', 'Sioux Falls SP','Sioux Falls Control','Custer SY','Custer SP','Custer Control'))

png("ITS_Pielou_genus.png", height = 4, width = 6, units = 'in', res = 600)
Pielou %>% 
  filter(Group.x %in% c('Sioux Falls SY', 'Sioux Falls SP', 'Sioux Falls Control', 'Custer SY', 'Custer SP', 'Custer Control')) %>%
  ggplot(aes(x=level_order, y=pielou, varwidth=TRUE, fill=factor(level_order))) + geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) + ggtitle("Alpha Diversity (Pielou's Evenness)") + labs(x="Group",y="Pielou Evenness") + theme_classic() + scale_color_manual(values=ITS_colors) + 
  scale_fill_manual(values=ITS_colors) + theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = 7, colour = "black"), axis.text.y = element_text(size = 7, colour = "black"))
dev.off()

shapiro.test(Pielou$pielou)
###p-value = 0.5857

##Run Anova based on Shapiro test
# Make sure Group is a factor
Pielou$Group.x <- factor(Pielou$Group.x)

# Run one-way ANOVA
anova_model_pielou <- aov(pielou ~ Group.x, data = Pielou)
summary(anova_model_pielou)
###########Df   Sum Sq   Mean Sq F value Pr(>F)  
#Group.x      5 0.008114 0.0016227   3.222 0.0174 *
#Residuals   34 0.017126 0.0005037 

TukeyHSD(anova_model_pielou)

# Residual normality
shapiro.test(residuals(anova_model_pielou))
##p-value = 0.457

# Visual check
qqnorm(residuals(anova_model_pielou)); qqline(residuals(anova_model_pielou), col = "red")

# Homogeneity of variance (Levene’s test)
install.packages("carData")
library(car)
leveneTest(pielou ~ Group.x, data = Pielou)
####Levene's Test for Homogeneity of Variance (center = median)
     ## Df F value Pr(>F)
###group  5  0.1309 0.9843
##34

Pielou$Group.x <- factor(Pielou$Group.x, levels = c(
  "Sioux Falls SY",
  "Sioux Falls SP",
  "Sioux Falls Control",
  "Custer SY",
  "Custer SP",
  "Custer Control"
))


###Plotting Pielou with significance
png("Pielou_genus_withstars.png", height = 5, width = 6, units = 'in', res = 600)
ggboxplot(Pielou, x = "Group.x", y = "pielou", fill = "Group.x", palette = ITS_colors) +
  stat_compare_means(method = "t.test", label = "p.signif", size = 4,
                     comparisons = list(
                       c("Sioux Falls SY", "Custer SY"),
                       c("Sioux Falls SP", "Custer SP"),
                       c("Sioux Falls Control", "Custer Control")
                     ),
                     label.y = c(0.785, 0.78)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(title = "Pielou's Evenness Index", y = "Pielou's Index", x = "Group")
dev.off()

