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

#############################################################################
##Repeat within drainage basins

colon_coords <- coords[1:57,]
chiriqui_coords<-coords[58:69,]

#Calculate ln(dij) using geodesic distances
col_geo_dist_matrix <- geodist(colon_coords, measure = "geodesic")  # in meters
col_log_geo_dist_matrix <- log(col_geo_dist_matrix + 1)
col_log_geo_dist <- as.dist(col_log_geo_dist_matrix)

ch_geo_dist_matrix <- geodist(chiriqui_coords, measure = "geodesic")  # in meters
ch_log_geo_dist_matrix <- log(ch_geo_dist_matrix + 1)
ch_log_geo_dist <- as.dist(ch_log_geo_dist_matrix)
#Calculate Euclidean genetic distances across all individual pairs
col_vcf <- read.vcfR("../marathrum_colon_filtered_final.recode.vcf")
col_genvcf <- vcfR2genlight(col_vcf)

ch_vcf <- read.vcfR("../marathrum_chiriqui_filtered_final.recode.vcf")
ch_genvcf <- vcfR2genlight(ch_vcf)

col_genetic_dist <- dist(col_genvcf, method = "euclidean")  # this gives a dist object
ch_genetic_dist <- dist(ch_genvcf, method = "euclidean")  # this gives a dist object

#Mantel test
col_mantel_result <- mantel(col_genetic_dist, col_log_geo_dist, method = "pearson", permutations = 999)
print(col_mantel_result)

ch_mantel_result <- mantel(ch_genetic_dist, ch_log_geo_dist, method = "pearson", permutations = 999)
print(ch_mantel_result)


#####################
####Plotting data####
#####################

gen_vec <- as.vector(ch_genetic_dist)
geo_vec <- as.vector(ch_log_geo_dist)

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

###################################
#####Repeat within each river######
###################################

#Geograohic distance matrixes
d_coords<-coords[1:26,]
a_coords<-coords[27:46,]
c_coords<-coords[47:57,]

d_geo_dist_matrix <- geodist(d_coords, measure = "geodesic")  # in meters
d_log_geo_dist_matrix <- log(d_geo_dist_matrix + 1)
d_log_geo_dist <- as.dist(d_log_geo_dist_matrix)

a_geo_dist_matrix <- geodist(a_coords, measure = "geodesic")  # in meters
a_log_geo_dist_matrix <- log(a_geo_dist_matrix + 1)
a_log_geo_dist <- as.dist(a_log_geo_dist_matrix)

c_geo_dist_matrix <- geodist(c_coords, measure = "geodesic")  # in meters
c_log_geo_dist_matrix <- log(c_geo_dist_matrix + 1)
c_log_geo_dist <- as.dist(c_log_geo_dist_matrix)


#Genetic distance matrixes
d_vcf <- read.vcfR("../marathrum_diego_filtered_final.recode.vcf")
d_genvcf <- vcfR2genlight(d_vcf)
d_genetic_dist <- dist(d_genvcf, method = "euclidean")  # this gives a dist object

a_vcf <- read.vcfR("../marathrum_aguacate_filtered_final.recode.vcf")
a_genvcf <- vcfR2genlight(a_vcf)
a_genetic_dist <- dist(a_genvcf, method = "euclidean")  # this gives a dist object

c_vcf <- read.vcfR("../marathrum_cocle_norte_filtered_final.recode.vcf")
c_genvcf <- vcfR2genlight(c_vcf)
c_genetic_dist <- dist(c_genvcf, method = "euclidean")  # this gives a dist object

#Mantel test
d_mantel_result <- mantel(d_genetic_dist, d_log_geo_dist, method = "pearson", permutations = 999)
print(d_mantel_result)

a_mantel_result <- mantel(a_genetic_dist, a_log_geo_dist, method = "pearson", permutations = 999)
print(a_mantel_result)

c_mantel_result <- mantel(c_genetic_dist, c_log_geo_dist, method = "pearson", permutations = 999)
print(c_mantel_result)

#####################
####Plotting data####
#####################

gen_vec <- as.vector(c_genetic_dist)
geo_vec <- as.vector(c_log_geo_dist)

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

