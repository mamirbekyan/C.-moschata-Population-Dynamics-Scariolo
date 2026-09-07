# Loading Libraries -----
library(vcfR)
library(poppr)
library(adegenet)
library(hierfstat)
library(ape)
library(Rcpp)
library(dartR)
library(adespatial)
library(ade4)
library(ggplot2)
library(PopGenReport)
library(plink)
library(tinytex)
library(RColorBrewer)
library(phytools)
library(radiator)
library(SNPRelate)
library(pegas)
library(TreeTools)
library(pheatmap)
library(vegan)
library(colorspace)
library(dplyr)
library(tibble)


# Setting Directory  ------

setwd("/home/bioinfo/Computing/stacks-2.65/research/Cucurbita/Cmax_rootstock/stacks1")
vcf <- read.vcfR("populations.haps.vcf", verbose = FALSE)
vcf_snp <- read.vcfR("populations.snps.vcf", verbose = FALSE)
pop <- read.table("/home/bioinfo/Computing/stacks-2.65/research/Cucurbita/Cmax_rootstock/popmaps/popmap1", header = FALSE, sep = "\t", stringsAsFactors = TRUE)

#loci <- vcfR2loci(vcf, return.alleles = TRUE)


#col <- c("red", "blue", "darkgreen", "purple")
gd <- vcfR::vcfR2genind(vcf, sep = "/", return.alleles = TRUE)
gd@pop <- as.factor(pop$V2)
gd@strata<- pop
gd_snp <- vcfR::vcfR2genind(vcf_snp, sep = "/", return.alleles = TRUE)
gd_snp@pop <- as.factor(pop$V2)
gd_snp@strata<- pop



# creating tree obj HAPS -----
gd_l <- missingno(gd, type = "loci", cutoff = 0.15, quiet=FALSE, freq=FALSE)
gd_g_l <- missingno(gd_l, type = "genotype", cutoff = 0.3, quiet=FALSE, freq=FALSE)

col <- divergingx_hcl(12, palette = "Spectral")

tree <- aboot(gd_g_l, tree = "upgma", distance = "nei.dist", sample = 1000, 
              showtree = TRUE, missing = "zero", mcutoff = 1, quiet = F, root = TRUE)
#tree_nj <- aboot(gd_g_l, tree = "nj", distance = "nei.dist", sample = 1000, 
#              showtree = TRUE, missing = "zero", mcutoff = 1, quiet = F, root = TRUE)


col <- c("#A71B4B", "#EB7803","blue" ,"#A08186" ,"orange", "darkgreen" ,"#584B9F")
plot.phylo(tree, type = "phylogram", cex = 0.9, font = 1, adj = 0, tip.color =  col[pop(gd_g_l)])
nodelabels(tree$node.label, adj = c(0, 0), frame = "n", cex = 0.6, font = 3, xpd = TRUE)
legend("topleft", legend = levels(gd_g_l$pop) , fill = col, col = col, border = "black", bty = "n", cex = 0.9, ncol = 1)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)
title(xlab = "Cucurbita maxima rootstock - Nei's Genetic Distance\nbased on 5409 haplotypic loci and 12569 variants")

# creating tree obj SNPS -----
gd_snp_l <- missingno(gd_snp, type = "loci", cutoff = 0.15, quiet=FALSE, freq=FALSE)
gd_snp_g_l <- missingno(gd_snp_l, type = "genotype", cutoff = 0.3, quiet=FALSE, freq=FALSE)


tree_snp <- aboot(gd_snp_g_l, tree = "upgma", distance = "nei.dist", sample = 1000, 
              showtree = TRUE, missing = "zero", mcutoff = 1, quiet = F, root = TRUE)

#col <- c("#A71B4B", "#D04939" ,"#EB7803", "#F5A736" ,"blue" ,"#A08186" ,"orange", "darkgreen" ,"#43CBB1", "palevioletred2" ,"#0080B2" ,"#584B9F")
plot.phylo(tree_snp, type = "phylogram", cex = 0.9, font = 1, adj = 0, tip.color =  col[pop(gd_snp_g_l)])
nodelabels(tree_snp$node.label, adj = c(0, 0), frame = "n", cex = 0.6, font = 3, xpd = TRUE)
legend("topleft", legend = levels(gd_snp_g_l$pop) , fill = col, col = col, border = "black", bty = "n", cex = 0.9, ncol = 1)
axisPhylo(side = 1, root.time = NULL, backward = TRUE)
title(xlab = "Cucurbita maxima rootstock - Nei's Genetic Distance\nbased on 20264 SNP loci and 40528 variants")





# (HAPS) Genetic Statistics (hierfstat) & WC_Fst comptuting and plotting -----
hier <- genind2hierfstat(gd_g_l, pop = gd_g_l$pop)


hier_samples <- genind2hierfstat(gd_g_l, pop = indNames(gd_g_l))

base_stat <- basic.stats(hier,diploid=TRUE,digits=3)
base_stat_samples <- basic.stats(hier_samples,diploid=TRUE,digits=3)

Fst <- pairwise.WCfst(hier, diploid = T)

Fst
cols_Fst <- divergingx_hcl(100, palette = "Spectral")
pheatmap(Fst, color = cols_Fst, cluster_rows=T, cluster_cols=T, na_col="white",main = "Pairwise Fst (HAPS)", display_numbers = T,
         number_format = "%.3f", number_color = "grey0", fontsize_number = 10)


# GD for Pops 
genepop_pop <- genind2genpop(gd_g_l, pop = gd_g_l@pop)
dist_pop <- dist.genpop(genepop_pop, method = 1, diag = T, upper = T)
dist_pop <- as.matrix(dist_pop)

cols_DistPop <- divergingx_hcl(100, palette = "Zissou 1")
pheatmap(dist_pop, color = cols_DistPop, cluster_rows = T,main = "Pairwise GD for Populations (HAPS)",
         cluster_cols = T, treeheight_row =0, treeheight_col =0, number_format = "%.3f", number_color = "black",  display_numbers = T)

# GD for Samples
genepop_sample <- genind2genpop(gd_g_l, pop = indNames(gd_g_l))
dist_sample <- dist.genpop(genepop_sample, method = 1, diag = T, upper = T)
dist_sample <- as.matrix(dist_sample)

cols_DistSamp <- divergingx_hcl(100, palette = "Temps")
pheatmap(dist_sample,color = cols_DistSamp, cluster_rows = T,clustering_distance_rows = "correlation", main = "Pairwise GD for Samples (HAPS)",
         cluster_cols = T, clustering_distance_cols = "correlation", treeheight_row =0, treeheight_col =0, number_format = "%.3f", number_color = "black",  display_numbers = F)


Ho_pop <- as.vector(base_stat["Ho"])
Hs_pop <- as.vector(base_stat["Hs"])
Fis_pop <- as.vector(base_stat["Fis"])
Overall_stats_pop <- t(data.frame(base_stat$overall))

base_stat_DF <- data.frame(Ho_pop, Hs_pop, Fis_pop)

Ho_DF <- as.data.frame(base_stat_samples[["Ho"]])


# (SNP) Genetic Statistics (hierfstat) & WC_Fst comptuting and plotting -----
hier_snp <- genind2hierfstat(gd_snp_g_l, pop = gd_snp_g_l$pop)


hier_samples_snp <- genind2hierfstat(gd_snp_g_l, pop = indNames(gd_snp_g_l))

base_stat_snp <- basic.stats(hier_snp,diploid=TRUE,digits=3)
base_stat_samples_snp <- basic.stats(hier_samples_snp,diploid=TRUE,digits=3)

Fst_snp <- pairwise.WCfst(hier_snp, diploid = T)

Fst_snp
cols_Fst <- divergingx_hcl(100, palette = "Spectral")
pheatmap(Fst_snp, color = cols_Fst, cluster_rows=T, cluster_cols=T, na_col="white",main = "Pairwise Fst (SNP)", display_numbers = T,
         number_format = "%.3f", number_color = "grey0", fontsize_number = 10)


# GD for Pops
genepop_pop_snp <- genind2genpop(gd_snp_g_l, pop = gd_snp_g_l@pop)
dist_pop_snp <- dist.genpop(genepop_pop_snp, method = 1, diag = T, upper = T)
#dist_pop_snp <- as.matrix(dist_pop)
dist_pop_snp <- as.matrix(dist_pop_snp)

pheatmap(dist_pop_snp, color = cols_DistPop, cluster_rows = T,main = "Pairwise GD for Populations (SNP)",
         cluster_cols = T, treeheight_row =0, treeheight_col =0, number_format = "%.3f", number_color = "black",  display_numbers = T)

# GD for Samples
genepop_sample_snp <- genind2genpop(gd_snp_g_l, pop = indNames(gd_snp_g_l))
dist_sample_snp <- dist.genpop(genepop_sample_snp, method = 1, diag = T, upper = T)
dist_sample_snp <- as.matrix(dist_sample_snp)

cols_DistSamp <- divergingx_hcl(100, palette = "Temps")
pheatmap(dist_sample_snp,color = cols_DistSamp, cluster_rows = T,clustering_distance_rows = "correlation", main = "Pairwise GD for Samples (SNP)",
         cluster_cols = T, clustering_distance_cols = "correlation", treeheight_row =0, treeheight_col =0, number_format = "%.3f", number_color = "black",  display_numbers = F)


Ho_pop_snp <- as.vector(base_stat_snp["Ho"])
Hs_pop_snp <- as.vector(base_stat_snp["Hs"])
Fis_pop_snp <- as.vector(base_stat_snp["Fis"])
Overall_stats_pop_snp <- t(data.frame(base_stat_snp$overall))

base_stat_DF_snp <- data.frame(Ho_pop_snp, Hs_pop_snp, Fis_pop_snp)

Ho_DF_snp <- as.data.frame(base_stat_samples_snp[["Ho"]])


### (HAPS) Create PCA ------

# Replace missing data with the mean allele frequencies
x = tab(gd_g_l, NA.method = "mean")
# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors haps(%)", ylim = c(0,25),
        names.arg = round(percent, 1))
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)
# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")
# Add a column containing individuals
ind_coords$Ind = indNames(gd_g_l)
# Add a column with the site IDs
ind_coords$Site = gd_g_l$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = col
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")
# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, linewidth=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15)
)
# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("PCA (HAPS)")+
  # custom theme
  ggtheme


#### (HAPS) Create DAPC -----

set.seed(5394)

crossval = xvalDapc(x, gd_g_l$pop, result = "groupMean", xval.plot = TRUE)


# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`

crossval$`Number of PCs Achieving Highest Mean Success`

crossval$`Number of PCs Achieving Lowest MSE`

numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)


# Run a DAPC using site IDs as priors
dapc1 = dapc(gd_g_l, gd_g_l$pop, n.pca = numPCs, n.da = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors hap(%)", ylim = c(0,60),
        names.arg = round(percent, 1))



#Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)
# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2")
# Add a column containing individuals
ind_coords$Ind = indNames(gd_g_l)
# Add a column with the site IDs
ind_coords$Site = gd_g_l$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = cols
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

#Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("DAPC HAP")+
  # custom theme
  ggtheme





### (SNP) Create PCA ------

# Replace missing data with the mean allele frequencies
x_snp = tab(gd_snp_g_l, NA.method = "mean")
# Perform PCA
pca1_snp = dudi.pca(x_snp, scannf = FALSE, scale = FALSE, nf = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent_snp = pca1_snp$eig/sum(pca1_snp$eig)*100
barplot(percent_snp, ylab = "Genetic variance explained by eigenvectors snp (%)", ylim = c(0,25),
        names.arg = round(percent_snp, 1))
# Create a data.frame containing individual coordinates
ind_coords_snp = as.data.frame(pca1_snp$li)
# Rename columns of dataframe
colnames(ind_coords_snp) = c("Axis1","Axis2","Axis3")
# Add a column containing individuals
ind_coords_snp$Ind = indNames(gd_snp_g_l)
# Add a column with the site IDs
ind_coords_snp$Site = gd_snp_g_l$pop
# Calculate centroid (average) position for each population
centroid_snp = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords_snp, FUN = mean)
# Add centroid coordinates to ind_coords dataframe
ind_coords_snp = left_join(ind_coords_snp, centroid_snp, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = col
# Custom x and y labels
xlab_snp = paste("Axis 1 (", format(round(percent_snp[1], 1), nsmall=1)," %)", sep="")
ylab_snp = paste("Axis 2 (", format(round(percent_snp[2], 1), nsmall=1)," %)", sep="")
# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, linewidth=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15)
)
# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords_snp, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid_snp, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab_snp, y = ylab_snp)+
  ggtitle("PCA (SNP)")+
  # custom theme
  ggtheme


#### (SNP) Create DAPC -----

set.seed(5394)

crossval_snp = xvalDapc(x_snp, gd_snp_g_l$pop, result = "groupMean", xval.plot = TRUE)


# Number of PCs with best stats (lower score = better)
crossval_snp$`Root Mean Squared Error by Number of PCs of PCA`

crossval_snp$`Number of PCs Achieving Highest Mean Success`

crossval_snp$`Number of PCs Achieving Lowest MSE`

numPCs_snp = as.numeric(crossval_snp$`Number of PCs Achieving Lowest MSE`)


# Run a DAPC using site IDs as priors
dapc1_snp = dapc(gd_snp_g_l, gd_snp_g_l$pop, n.pca = numPCs_snp, n.da = 3)
# Analyse how much percent of genetic variance is explained by each axis
percent_snp = dapc1_snp$eig/sum(dapc1_snp$eig)*100
barplot(percent_snp, ylab = "Genetic variance explained by eigenvectors (%) snp", ylim = c(0,60),
        names.arg = round(percent, 1))



#Create a data.frame containing individual coordinates
ind_coords_snp = as.data.frame(dapc1_snp$ind.coord)
# Rename columns of dataframe
colnames(ind_coords_snp) = c("Axis1","Axis2")
# Add a column containing individuals
ind_coords_snp$Ind = indNames(gd_snp_g_l)
# Add a column with the site IDs
ind_coords_snp$Site = gd_snp_g_l$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords_snp, FUN = mean)
# Add centroid coordinates to ind_coords_snp dataframe
ind_coords_snp = left_join(ind_coords_snp, centroid, by = "Site", suffix = c("",".cen"))
# Define colour palette
cols = cols
# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

#Scatter plot axis 1 vs. 2
ggplot(data = ind_coords_snp, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("DAPC SNP")+
  # custom theme
  ggtheme


## (HAPS) Export_tables and file for other analyses (STRUCTURE) -----

setwd("/home/bioinfo/Computing/stacks-2.65/research/Cucurbita/Cmax_rootstock/stacks1/Genetic_Statistics")

# Tables
write.csv2(Ho_DF, file = "Ho_samples.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(dist_sample, file = "Nei_dist_samples.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(dist_pop , file = "Nei_dist_pop.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(Fst , file = "Fst_pairwise.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )

write.csv2(base_stat_DF , file = "base_stat_DF.csv", append = T, quote = F, sep = ",",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(Overall_stats_pop , file = "Overall_stats_pop.csv", append = T, quote = F, sep = ",",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )

write.csv2(Ho_pop, file = "Ho.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )


## (SNP) Export_tables and file for other analyses (STRUCTURE) -----



# Tables
write.csv2(Ho_DF_snp, file = "Ho_samples_snp.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(dist_sample_snp, file = "Nei_dist_samples_snp.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(dist_pop_snp , file = "Nei_dist_pop_snp.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(Fst_snp , file = "Fst_pairwise_snp.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )

write.csv2(base_stat_DF_snp , file = "base_stat_DF_snp.csv", append = T, quote = F, sep = ",",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )
write.csv2(Overall_stats_pop_snp , file = "Overall_stats_pop_snp.csv", append = T, quote = F, sep = ",",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )



write.csv2(Ho_pop_snp, file = "Ho_pop_snp.csv", append = FALSE, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )




#iQTree input creation
loci <- vcfR2loci(vcf,  return.alleles = F)

tab <- gd_g_l@tab
tab <- t(tab)
write.csv2(tab , file = "loci.csv", append = T, quote = F, sep = ";",
           eol = "\n", na = "NA", dec = ".", row.names = TRUE,
           col.names = TRUE, qmethod = c("escape", "double"),
           fileEncoding = "" )



# AMOVA computing ----

Ind <- as.vector(indNames(gd_g_l))
Pop <- as.vector(gd_g_l@pop)
Species <- as.vector(pop$V4)
pop1 <- data.frame(Ind, Pop, Species)

gd_g_l@strata <- pop1



amova <- poppr.amova(gd_g_l, ~V3,  clonecorrect = F,
                     within = T,
                     dist = NULL,
                     squared = F,
                     freq = T,
                     correction = "quasieuclid",
                     filter = FALSE,
                     threshold = 0,
                     algorithm = "average_neighbor",
                     threads = 1L,
                     missing = "loci",
                     cutoff = 1,
                     quiet = FALSE,
                     method = "ade4",
                     nperm = 999)

print(amova)

# Test for significance
set.seed(1234)
amova.test <- randtest(amova, nrepet = 999) 
print(amova)
print(amova.test)
plot(amova.test)



### STRUCTURE input creation -------------
write.struct(hier_samples,ilab=Ind,pop=NULL,MARKERNAMES=T,MISSING=-9,fname="fennel_filtered.str")



