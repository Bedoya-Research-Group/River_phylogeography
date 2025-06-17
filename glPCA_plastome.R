library(adegenet)
library(vcfR)
library(ggplot2)

setwd("/Users/abedoya/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/")

vcf <- read.vcfR("Marathrum_samples_filtered_Panama.vcf.recode.vcf")
#Read as genlight object
genlight_obj <- vcfR2genlight(vcf)

#check indiv names
vcf@gt
colnames(vcf@gt)
current_names <- colnames(vcf@gt)[-1]
print(current_names)

new_names <- sub("_plastome.*", "", current_names)
colnames(vcf@gt)[-1] <- new_names

pop_vector <- c(
  rep("D", 26),
  rep("A", 20),
  rep("C", 11),
  rep("Pa", 10),
  rep("Pi", 2)
)


#Assign pop vector to genlight object and check
pop(genlight_obj) <- pop_vector
head(pop(genlight_obj))

pca_result <- glPca(genlight_obj, nf=3)

# Get PC scores
pca_scores <- as.data.frame(pca_result$scores)

# Add population info
pca_scores$pop <- pop(genlight_obj)

# Extract indiv. names
ind_names <- indNames(genlight_obj)
pca_scores$ind <- ind_names
#pca_scores$ind <- ind_names
pca_scores$label <- ifelse(pca_scores$pop == "A", pca_scores$ind, NA)


ggplot(pca_scores, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 3) +
  labs(title = "PCA of SNP data", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank())

pca_table <- pca_scores[, c("ind", "pop", "PC1", "PC2", "PC3")]

pca_table

##################
#######COCLE######
##################

vcf2 <- read.vcfR("Marathrum_samples_filtered_Cocle.vcf.recode.vcf")
#Read as genlight object
genlight_obj <- vcfR2genlight(vcf2)

#check indiv names
vcf2@gt
colnames(vcf2@gt)
current_names <- colnames(vcf2@gt)[-1]
print(current_names)

new_names <- sub("_plastome.*", "", current_names)
colnames(vcf@gt)[-1] <- new_names

pop_vector <- c(
  rep("D", 26),
  rep("A", 20),
  rep("C", 11)
)


#Assign pop vector to genlight object and check
pop(genlight_obj) <- pop_vector
head(pop(genlight_obj))

pca_result <- glPca(genlight_obj, nf=3)

# Get PC scores
pca_scores <- as.data.frame(pca_result$scores)

# Add population info
pca_scores$pop <- pop(genlight_obj)

# Extract indiv. names
ind_names <- indNames(genlight_obj)
pca_scores$ind <- ind_names
#pca_scores$ind <- ind_names
pca_scores$label <- ifelse(pca_scores$pop == "A", pca_scores$ind, NA)


ggplot(pca_scores, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 3) +
  labs(title = "PCA of SNP data", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank())

pca_table <- pca_scores[, c("ind", "pop", "PC1", "PC2", "PC3")]

pca_table
