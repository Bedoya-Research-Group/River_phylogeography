##Gene flow and genetic diversity across space

library(vegan)
library(geodist)
library(vcfR)
library(adegenet)
library(ggplot2)

setwd("~setwd("~setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/River_phylogeography/")
coords <- read.csv("coordinates.csv", row.names = 1)

#Calculate ln(dij) using geodesic distances
geo_dist_matrix <- geodist(coords, measure = "geodesic")  # in meters
log_geo_dist_matrix <- log(geo_dist_matrix + 1)
log_geo_dist <- as.dist(log_geo_dist_matrix)

#Calculate Euclidean genetic distances across all individual pairs
vcf <- read.vcfR("../marathrum_panama_filtered_final.recode.vcf")
genvcf <- vcfR2genlight(vcf)

genetic_dist <- dist(genvcf, method = "euclidean")  # this gives a dist object


#Mantel test
mantel_result <- mantel(genetic_dist, log_geo_dist, method = "pearson", permutations = 999)
print(mantel_result)

#####################
####Plotting data####
#####################

gen_vec <- as.vector(genetic_dist)
geo_vec <- as.vector(log_geo_dist)

# Combine into a data frame
dist_df <- data.frame(
  Genetic = gen_vec,
  Log_Geographic = geo_vec
)


ggplot(dist_df, aes(x = Log_Geographic, y = Genetic)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  theme_minimal() +
  labs(
    x = "ln(dij)",
    y = "Genetic Distance"
  )
cor(dist_df$Genetic, dist_df$Log_Geographic, method = "pearson")



