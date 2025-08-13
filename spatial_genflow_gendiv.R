##Gene flow and genetic diversity across space

library(vegan)
library(geodist)
library(vcfR)
library(adegenet)
library(ggplot2)
library(hierfstat)
library(ade4)
library(dplyr)

setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/River_phylogeography/")
coords <- read.csv("coordinates.csv", row.names = 1)

#Calculate ln(dij) using geodesic distances
geo_dist_matrix <- geodist(coords, measure = "geodesic")  # in meters
log_geo_dist_matrix <- log(geo_dist_matrix + 1)
log_geo_dist <- as.dist(log_geo_dist_matrix)

#Calculate Euclidean genetic distances across all individual pairs
vcf <- read.vcfR("../marathrum_panama_filtered_final.recode.vcf")
#vcf <- read.vcfR("Data_plastome/Marathrum_samples_filtered_Panama.vcf.recode.vcf") for plastome analysis
genvcf <- vcfR2genlight(vcf)

genetic_dist <- dist(genvcf, method = "euclidean")


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

individual_ids<-read.csv("colon_sample_ids.txt",head=F)
rownames(colon_coords) <- individual_ids[1:57,]

#Calculate ln(dij) using geodesic distances
col_geo_dist_matrix <- geodist(colon_coords, measure = "geodesic")  # in meters
col_log_geo_dist_matrix <- log(col_geo_dist_matrix + 1)
col_log_geo_dist <- as.dist(col_log_geo_dist_matrix)

ch_geo_dist_matrix <- geodist(chiriqui_coords, measure = "geodesic")  # in meters
ch_log_geo_dist_matrix <- log(ch_geo_dist_matrix + 1)
ch_log_geo_dist <- as.dist(ch_log_geo_dist_matrix)
#Calculate Euclidean genetic distances across all individual pairs
col_vcf <- read.vcfR("../marathrum_colon_filtered_final.recode.vcf")
#col_vcf <- read.vcfR("Data_plastome/Marathrum_samples_filtered_Cocle.vcf.recode.vcf")#for plastome data
col_genvcf <- vcfR2genlight(col_vcf)


ch_vcf <- read.vcfR("../marathrum_chiriqui_filtered_final.recode.vcf")
#ch_vcf <- read.vcfR("Data_plastome/marathrum_chiriqui_filtered_final.recode.vcf")#for plastome analysis
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

#Geographic distance matrices
d_coords<-coords[1:26,]
a_coords<-coords[27:46,]
c_coords<-coords[47:57,]
pa_coords<-coords[58:67,]

d_geo_dist_matrix <- geodist(d_coords, measure = "geodesic")  # in meters
d_log_geo_dist_matrix <- log(d_geo_dist_matrix + 1)
d_log_geo_dist <- as.dist(d_log_geo_dist_matrix)

a_geo_dist_matrix <- geodist(a_coords, measure = "geodesic")  # in meters
a_log_geo_dist_matrix <- log(a_geo_dist_matrix + 1)
a_log_geo_dist <- as.dist(a_log_geo_dist_matrix)

c_geo_dist_matrix <- geodist(c_coords, measure = "geodesic")  # in meters
c_log_geo_dist_matrix <- log(c_geo_dist_matrix + 1)
c_log_geo_dist <- as.dist(c_log_geo_dist_matrix)

pa_geo_dist_matrix <- geodist(pa_coords, measure = "geodesic")  # in meters
pa_log_geo_dist_matrix <- log(pa_geo_dist_matrix + 1)
pa_log_geo_dist <- as.dist(pa_log_geo_dist_matrix)
############################################################
##IF USING DISTANCES ESTIMATED FROM RIVER TRAJECTORY
# Convert geodist objects to matrices, then to data frames
#c_geo_dist_df <- as.data.frame(as.matrix(c_geo_dist_matrix))
# Save to CSV
write.csv(c_geo_dist_df, file = "geo_dist_matrix_c.csv", quote = FALSE)
#Modify to convert matrices into distances estimated from river trajectory per river
# read matrices

a_path <- read.csv("geo_dist_matrix_a.csv", row.names = 1)
a_path_dist_matrix <- as.matrix(a_path)
a_log_path_dist_matrix <- log(a_path_dist_matrix + 1)
a_log_path_dist <- as.dist(a_log_path_dist_matrix)

d_path <- read.csv("geo_dist_matrix_d.csv", row.names = 1)
d_path_dist_matrix <- as.matrix(d_path)
d_log_path_dist_matrix <- log(d_path_dist_matrix + 1)
d_log_path_dist <- as.dist(d_log_path_dist_matrix)

c_path <- read.csv("geo_dist_matrix_c.csv", row.names = 1)
c_path_dist_matrix <- as.matrix(c_path)
c_log_path_dist_matrix <- log(c_path_dist_matrix + 1)
c_log_path_dist <- as.dist(c_log_path_dist_matrix)


#Genetic distance matrices
d_vcf <- read.vcfR("../marathrum_diego_filtered_final.recode.vcf")
#d_vcf <- read.vcfR("Data_plastome/marathrum_diego_filtered_final.recode.vcf")#for plastom analysis
d_genvcf <- vcfR2genlight(d_vcf)
d_genetic_dist <- dist(d_genvcf, method = "euclidean")  # this gives a dist object

a_vcf <- read.vcfR("../marathrum_aguacate_filtered_final.recode.vcf")
#a_vcf <- read.vcfR("Data_plastome/marathrum_aguacate_filtered_final.recode.vcf") #for plastome analysis
a_genvcf <- vcfR2genlight(a_vcf)
a_genetic_dist <- dist(a_genvcf, method = "euclidean")  # this gives a dist object

c_vcf <- read.vcfR("../marathrum_cocle_norte_filtered_final.recode.vcf")
#c_vcf <- read.vcfR("Data_plastome/marathrum_cocle_norte_filtered_final.recode.vcf")#for plastome analysis
c_genvcf <- vcfR2genlight(c_vcf)
c_genetic_dist <- dist(c_genvcf, method = "euclidean")  # this gives a dist object

pa_vcf <- read.vcfR("../marathrum_paraiso_filtered_final.recode.vcf")
#pa_vcf <- read.vcfR("Data_plastome/marathrum_paraiso_filtered_final.recode.vcf")#for plastome analysis
pa_genvcf <- vcfR2genlight(pa_vcf)
pa_genetic_dist <- dist(pa_genvcf, method = "euclidean")  # this gives a dist object

#Mantel test
d_mantel_result <- mantel(d_genetic_dist, d_log_geo_dist, method = "pearson", permutations = 999)
print(d_mantel_result)

a_mantel_result <- mantel(a_genetic_dist, a_log_geo_dist, method = "pearson", permutations = 999)
print(a_mantel_result)

c_mantel_result <- mantel(c_genetic_dist, c_log_geo_dist, method = "pearson", permutations = 999)
print(c_mantel_result)

pa_mantel_result <- mantel(pa_genetic_dist, pa_log_geo_dist, method = "pearson", permutations = 999)
print(pa_mantel_result)

d_mantel_result_path <- mantel(d_genetic_dist, d_log_path_dist, method = "pearson", permutations = 999)
print(d_mantel_result)

a_mantel_result_path <- mantel(a_genetic_dist, a_log_path_dist, method = "pearson", permutations = 999)
print(a_mantel_result)

c_mantel_result_path <- mantel(c_genetic_dist, c_log_path_dist, method = "pearson", permutations = 999)
print(c_mantel_result)

#Cannot do for Paraiso because of damn building and lots of modification to the land that has impacted the course of the river
#pa_mantel_result_path <- mantel(pa_genetic_dist, pa_log_path_dist, method = "pearson", permutations = 999)
#print(pa_mantel_result)

#####################
####Plotting data####
#####################

gen_vec <- as.vector(pa_genetic_dist)
geo_vec <- as.vector(pa_log_geo_dist)

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



####################################################################
###To test if genetic distance is stronger within or across rivers##
########(as a proxy for dispersal within or across rivers)##########
####################################################################

#all
#Calculate ln(dij) using geodesic distances

individual_ids_panama<-read.csv("panama_sample_ids.txt",head=F)
rownames(coords) <- individual_ids_panama[,]
geo_dist_matrix <- geodist(coords, measure = "geodesic")  # in meters
log_geo_dist_matrix <- log(geo_dist_matrix + 1)
rownames(log_geo_dist_matrix) <- rownames(coords)
colnames(log_geo_dist_matrix) <- rownames(coords)
log_geo_dist <- as.dist(log_geo_dist_matrix)

#To add same labels to the genetic matrix:
individual_ids_panama <- individual_ids_panama$V1
pa_genetic_mat <- as.matrix(genetic_dist)
rownames(pa_genetic_mat) <- individual_ids_panama
colnames(pa_genetic_mat) <- individual_ids_panama
pa_genetic_dist <- as.dist(pa_genetic_mat)

#names(river_labels) <- individual_ids

individuals <- labels(pa_genetic_dist)

river <- ifelse(grepl("Diego", individuals), "Diego",
                ifelse(grepl("Aguacate", individuals), "Aguacate",
                       ifelse(grepl("Cocle", individuals), "Cocle",
                              ifelse(grepl("Paraiso", individuals), "Paraiso",NA))))

stopifnot(all(labels(pa_genetic_dist) == labels(log_geo_dist))) #checking order in geo and gen distances matrices is the same

###########
#Mantel test across all panama(repeat from above but here for clarity)
###########

#mantel_overall <- mantel(pa_genetic_dist, log_geo_dist, method = "pearson", permutations = 999)
#print(mantel_overall)

# Function to extract and test per river
#mantel_within_river <- function(river_name) {
#  inds <- individuals[river == river_name]
#  gen_sub <- as.dist(as.matrix(pa_genetic_dist)[inds, inds])
#  geo_sub <- as.dist(as.matrix(log_geo_dist)[inds, inds])
#  mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
#}

#mantel_diego <- mantel_within_river("Diego")
#mantel_aguacate <- mantel_within_river("Aguacate")
#mantel_cocle <- mantel_within_river("Cocle")
#mantel_paraiso <- mantel_within_river("Paraiso")

#print(mantel_diego)
#print(mantel_aguacate)
#print(mantel_cocle)
#print(mantel_paraiso)

# River identity dissimilarity matrix: 0 = same river, 1 = different
river_matrix <- outer(river, river, FUN = function(x, y) as.numeric(x != y))
river_dist <- as.dist(river_matrix)

# Partial Mantel test
partial_mantel <- mantel.partial(pa_genetic_dist, log_geo_dist, river_dist, method = "pearson", permutations = 999)
print(partial_mantel)

##Plotting
# Create pairwise comparison data frame
plot_df <- data.frame(
  GeneticDist = as.vector(as.matrix(pa_genetic_dist)),
  GeoDist = as.vector(as.matrix(log_geo_dist)),
  River1 = rep(river, each = length(river)),
  River2 = rep(river, times = length(river))
)

# Label comparisons
plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, "Within", "Between")

# Plot: color by Within vs. Between
ggplot(plot_df, aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Within" = "grey", "Between" = "black")) +
  labs(x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal() +
  theme(legend.title = element_blank())
#########################
##Colon
rownames(col_log_geo_dist_matrix) <- rownames(colon_coords)
colnames(col_log_geo_dist_matrix) <- rownames(colon_coords)

col_log_geo_dist <- as.dist(col_log_geo_dist_matrix)

#To add same labels to the genetic matrix:
individual_ids <- individual_ids$V1
col_genetic_mat <- as.matrix(col_genetic_dist)
rownames(col_genetic_mat) <- individual_ids
colnames(col_genetic_mat) <- individual_ids
col_genetic_dist <- as.dist(col_genetic_mat)

#names(river_labels) <- individual_ids

individuals <- labels(col_genetic_dist)

river <- ifelse(grepl("Diego", individuals), "Diego",
                ifelse(grepl("Aguacate", individuals), "Aguacate",
                       ifelse(grepl("Cocle", individuals), "Cocle", NA)))

stopifnot(all(labels(col_genetic_dist) == labels(col_log_geo_dist))) #checking order in geo and gen distances matrices is the same

#Mantel test across all colon(repeat from above but here for clarity)
mantel_overall <- mantel(col_genetic_dist, col_log_geo_dist, method = "pearson", permutations = 999)

# Function to extract and test per river
mantel_within_river <- function(river_name) {
  inds <- individuals[river == river_name]
  gen_sub <- as.dist(as.matrix(col_genetic_dist)[inds, inds])
  geo_sub <- as.dist(as.matrix(col_log_geo_dist)[inds, inds])
  mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
}

#mantel_diego <- mantel_within_river("Diego")
#mantel_aguacate <- mantel_within_river("Aguacate")
#mantel_cocle <- mantel_within_river("Cocle")

#print(mantel_diego)
#print(mantel_aguacate)
#print(mantel_cocle)

# River identity dissimilarity matrix: 0 = same river, 1 = different
river_matrix <- outer(river, river, FUN = function(x, y) as.numeric(x != y))
river_dist <- as.dist(river_matrix)

# Partial Mantel test
partial_mantel <- mantel.partial(col_genetic_dist, col_log_geo_dist, river_dist, method = "pearson", permutations = 999)
print(partial_mantel)

##Plotting
# Create pairwise comparison data frame
plot_df <- data.frame(
  GeneticDist = as.vector(as.matrix(col_genetic_dist)),
  GeoDist = as.vector(as.matrix(col_log_geo_dist)),
  River1 = rep(river, each = length(river)),
  River2 = rep(river, times = length(river))
)

# Label comparisons
plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, "Within", "Between")

# Plot: color by Within vs. Between
ggplot(plot_df, aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Within" = "grey", "Between" = "black")) +
  labs(x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal() +
  theme(legend.title = element_blank())

##Chiriqui
rownames(ch_log_geo_dist_matrix) <- rownames(chiriqui_coords)
colnames(ch_log_geo_dist_matrix) <- rownames(chiriqui_coords)

ch_log_geo_dist <- as.dist(ch_log_geo_dist_matrix)
ch_mat <- as.matrix(ch_log_geo_dist)

current_names <- rownames(ch_mat) #current names
# Modify the last two
current_names[(length(current_names)-1):length(current_names)] <- c("Piedra_01", "Piedra_02")
# Apply new names
rownames(ch_mat) <- colnames(ch_mat) <- current_names

#To add same labels to the genetic matrix:
individual_ids_chiriqui <- c("Paraiso_01","Paraiso_02","Paraiso_03","Paraiso_04",
                             "Paraiso_05","Paraiso_06","Paraiso_07","Paraiso_08",
                             "Paraiso_09","Paraiso_10","Piedra_01","Piedra_02")
ch_genetic_mat <- as.matrix(ch_genetic_dist)
rownames(ch_genetic_mat) <- individual_ids_chiriqui
colnames(ch_genetic_mat) <- individual_ids_chiriqui
ch_genetic_dist <- as.dist(ch_genetic_mat)

ch_log_geo_dist <- as.dist(ch_mat) #Convert back to dist object
#names(river_labels) <- individual_ids

individuals <- labels(ch_genetic_dist)

river <- ifelse(grepl("Paraiso", individuals), "Paraiso",
                ifelse(grepl("Piedra", individuals), "Piedra", NA))

stopifnot(all(labels(ch_genetic_dist) == labels(ch_log_geo_dist))) #checking order in geo and gen distances matrices is the same

# Function to extract and test per river
mantel_within_river <- function(river_name) {
  inds <- individuals[river == river_name]
  gen_sub <- as.dist(as.matrix(ch_genetic_dist)[inds, inds])
  geo_sub <- as.dist(as.matrix(ch_log_geo_dist)[inds, inds])
  mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
}

#mantel_diego <- mantel_within_river("Diego")
#mantel_aguacate <- mantel_within_river("Aguacate")
#mantel_cocle <- mantel_within_river("Cocle")

#print(mantel_diego)
#print(mantel_aguacate)
#print(mantel_cocle)

# River identity dissimilarity matrix: 0 = same river, 1 = different
river_matrix <- outer(river, river, FUN = function(x, y) as.numeric(x != y))
river_dist <- as.dist(river_matrix)

# Partial Mantel test
partial_mantel <- mantel.partial(ch_genetic_dist, ch_log_geo_dist, river_dist, method = "pearson", permutations = 999)
print(partial_mantel)

##Plotting
# Create pairwise comparison data frame
plot_df <- data.frame(
  GeneticDist = as.vector(as.matrix(ch_genetic_dist)),
  GeoDist = as.vector(as.matrix(ch_log_geo_dist)),
  River1 = rep(river, each = length(river)),
  River2 = rep(river, times = length(river))
)

# Label comparisons
plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, "Within", "Between")

# Plot: color by Within vs. Between
ggplot(plot_df, aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Within" = "grey", "Between" = "black")) +
  labs(x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal() +
  theme(legend.title = element_blank())


#plot_df <- data.frame(
#  GeneticDist = as.vector(as.matrix(col_genetic_dist)),
#  GeoDist = as.vector(as.matrix(col_log_geo_dist)),
#  River1 = rep(river, each = length(river)),
#  River2 = rep(river, times = length(river))
#)

#plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, plot_df$River1, "Between")

#ggplot(plot_df[plot_df$Comparison != "Between",], aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
#  geom_point(alpha = 0.5) +
#  geom_smooth(method = "lm") +
#  theme_minimal()

####################################################
######Doing the same between adjacent rivers########
####################################################
#####Distances between adjacent rivers######
####################################################
ad_coords<-coords[1:46,]
ac_coords<-coords[27:57,]

#Geographic distances between A-D
ad_geo_dist_matrix <- geodist(ad_coords, measure = "geodesic")  # in meters
ad_log_geo_dist_matrix <- log(ad_geo_dist_matrix + 1)
ad_log_geo_dist <- as.dist(ad_log_geo_dist_matrix)

#Geographic distances between A-C
ac_geo_dist_matrix <- geodist(ac_coords, measure = "geodesic")  # in meters
ac_log_geo_dist_matrix <- log(ac_geo_dist_matrix + 1)
ac_log_geo_dist <- as.dist(ac_log_geo_dist_matrix)

##Genetic Distances
ad_vcf <- read.vcfR("../marathrum_aguacate_diego_filtered_final.recode.vcf")
ad_genvcf <- vcfR2genlight(ad_vcf)
ad_genetic_dist <- dist(ad_genvcf, method = "euclidean")  # this gives a dist object

ac_vcf <- read.vcfR("../marathrum_aguacate_cocle_filtered_final.recode.vcf")
ac_genvcf <- vcfR2genlight(ac_vcf)
ac_genetic_dist <- dist(ac_genvcf, method = "euclidean")  # this gives a dist object

##Steps for partial mantel
  ##AD
rownames(ad_log_geo_dist_matrix) <- rownames(colon_coords[1:46,])
colnames(ad_log_geo_dist_matrix) <- rownames(colon_coords[1:46,])

ad_log_geo_dist <- as.dist(ad_log_geo_dist_matrix)

  ##AC
rownames(ac_log_geo_dist_matrix) <- rownames(colon_coords[27:57,])
colnames(ac_log_geo_dist_matrix) <- rownames(colon_coords[27:57,])

ac_log_geo_dist <- as.dist(ac_log_geo_dist_matrix)

#To add same labels to the genetic matrix:
  
  ##AD
ad_individual_ids <- individual_ids[1:46]
ad_genetic_mat <- as.matrix(ad_genetic_dist)
rownames(ad_genetic_mat) <- ad_individual_ids
colnames(ad_genetic_mat) <- ad_individual_ids
ad_genetic_dist <- as.dist(ad_genetic_mat)

#names(river_labels) <- individual_ids
individuals <- labels(ad_genetic_dist)

river <- ifelse(grepl("Diego", individuals), "Diego",
                ifelse(grepl("Aguacate", individuals), "Aguacate", NA))

stopifnot(all(labels(ad_genetic_dist) == labels(ad_log_geo_dist))) #checking order in geo and gen distances matrices is the same

  ##AC
ac_individual_ids <- individual_ids[27:57]
ac_genetic_mat <- as.matrix(ac_genetic_dist)
rownames(ac_genetic_mat) <- ac_individual_ids
colnames(ac_genetic_mat) <- ac_individual_ids
ac_genetic_dist <- as.dist(ac_genetic_mat)

###########################################
  ##AD
individuals <- labels(ad_genetic_dist)

river <- ifelse(grepl("Diego", individuals), "Diego",
                ifelse(grepl("Aguacate", individuals), "Aguacate", NA))

stopifnot(all(labels(ad_genetic_dist) == labels(ad_log_geo_dist))) #checking order in geo and gen distances matrices is the same

# Function to extract and test per river
mantel_within_river <- function(river_name) {
  inds <- individuals[river == river_name]
  gen_sub <- as.dist(as.matrix(ad_genetic_dist)[inds, inds])
  geo_sub <- as.dist(as.matrix(ad_log_geo_dist)[inds, inds])
  mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
}

# River identity dissimilarity matrix: 0 = same river, 1 = different
river_matrix <- outer(river, river, FUN = function(x, y) as.numeric(x != y))
river_dist <- as.dist(river_matrix)
###########################################
  ##AC
individuals <- labels(ac_genetic_dist)

river <- ifelse(grepl("Aguacate", individuals), "Aguacate",
                ifelse(grepl("Cocle", individuals), "Cocle", NA))

stopifnot(all(labels(ac_genetic_dist) == labels(ac_log_geo_dist))) #checking order in geo and gen distances matrices is the same

# Function to extract and test per river
mantel_within_river <- function(river_name) {
  inds <- individuals[river == river_name]
  gen_sub <- as.dist(as.matrix(ac_genetic_dist)[inds, inds])
  geo_sub <- as.dist(as.matrix(ac_log_geo_dist)[inds, inds])
  mantel(gen_sub, geo_sub, method = "pearson", permutations = 999)
}

# River identity dissimilarity matrix: 0 = same river, 1 = different
river_matrix <- outer(river, river, FUN = function(x, y) as.numeric(x != y))
river_dist <- as.dist(river_matrix)

###########################################
# Partial Mantel test
  ##AD
partial_mantel <- mantel.partial(ad_genetic_dist, ad_log_geo_dist, river_dist, method = "pearson", permutations = 999)
print(partial_mantel)

  ##AC
partial_mantel <- mantel.partial(ac_genetic_dist, ac_log_geo_dist, river_dist, method = "pearson", permutations = 999)
print(partial_mantel)

##Plotting

########################
  ##AD
# Create pairwise comparison data frame
plot_df <- data.frame(
  GeneticDist = as.vector(as.matrix(ad_genetic_dist)),
  GeoDist = as.vector(as.matrix(ad_log_geo_dist)),
  River1 = rep(river, each = length(river)),
  River2 = rep(river, times = length(river))
)

# Label comparisons
plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, "Within", "Between")

# Plot: color by Within vs. Between
ggplot(plot_df, aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Within" = "grey", "Between" = "black")) +
  labs(x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal() +
  theme(legend.title = element_blank())

########################
  ##AC
# Create pairwise comparison data frame
plot_df <- data.frame(
  GeneticDist = as.vector(as.matrix(ac_genetic_dist)),
  GeoDist = as.vector(as.matrix(ac_log_geo_dist)),
  River1 = rep(river, each = length(river)),
  River2 = rep(river, times = length(river))
)

# Label comparisons
plot_df$Comparison <- ifelse(plot_df$River1 == plot_df$River2, "Within", "Between")

# Plot: color by Within vs. Between
ggplot(plot_df, aes(x = GeoDist, y = GeneticDist, color = Comparison)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Within" = "grey", "Between" = "black")) +
  labs(x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal() +
  theme(legend.title = element_blank())

######################
#Plotting only the "between" for clarity
plot_df_between <- subset(plot_df, Comparison == "Between")

# Plot only the between-river points
ggplot(plot_df_between, aes(x = GeoDist, y = GeneticDist)) +
  geom_point(alpha = 0.6, color = "black") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = "Genetic vs Geographic Distance (Between Rivers Only)",
       x = "Geographic Distance (log meters)",
       y = "Genetic Distance (Euclidean)") +
  theme_minimal()

#########################################################
#######  Diversity Stats within rivers #######
#########################################################
setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/")
#genind_obj <- vcfR2genind(col_vcf)
vcf_stats <- read.vcfR("marathrum4diversity_stats.recode.vcf")
genind_obj <- vcfR2genind(vcf_stats)
#Populations assignment

pop(genind_obj) <- factor(c(
  rep("Diego_down", 9),
  rep("Diego_mid", 7),
  rep("Diego_up", 10),
  rep("Aguacate_down", 8),
  rep("Aguacate_mid", 5),
  rep("Aguacate_up", 7),
  rep("Paraiso_up", 5),
  rep("Paraiso_down", 5)
))

geno <- tab(genind_obj, NA.method = "mean")

##Calculating Ho with adegenet

X <- tab(genind_obj, freq = FALSE, NA.method = "asis") #ignoring missing sites
loc_factor <- rep(locNames(genind_obj), times = nAll(genind_obj))

calc_ind_Ho <- function(ind_row) {
  split_loci <- split(ind_row, loc_factor)
  het_vec <- vapply(split_loci, function(vec) {
    if (any(is.na(vec))) return(NA_real_)
    nonzero <- vec[vec > 0]
    if (length(nonzero) >= 2) 1 else 0
  }, numeric(1))
  mean(het_vec, na.rm = TRUE)
}

ind_Ho <- apply(X, 1, calc_ind_Ho)

# Combine with population info
df_indHo <- data.frame(
  Individual = indNames(genind_obj),
  Population = pop(genind_obj),
  Ho = ind_Ho
)

df_indHo %>%
  group_by(Population) %>%
  summarise(
    Mean_Ho = mean(Ho, na.rm = TRUE),
    SD_Ho = sd(Ho, na.rm = TRUE),
    N = n()
  )
#Plot
df_indHo$Population <- factor( #reordering populations
  df_indHo$Population,
  levels = c(
    "Diego_up", "Diego_mid", "Diego_down",
    "Aguacate_up", "Aguacate_mid", "Aguacate_down",
    "Paraiso_up", "Paraiso_down"
  )
)

ggplot(df_indHo, aes(x = Population, y = Ho, fill = Population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  theme_minimal() +
  labs(
    y = "Observed Heterozygosity (Ho)",
    x = NULL,
    title = "Observed Heterozygosity per Population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )



# Function to calculate Ho for one individual (wiht imputation so the method above is better)
#calculate_Ho <- function(row) {
#  het_count <- 0
#  for (i in seq(1, ncol(geno), by = 2)) {
#    alleles <- row[i:(i + 1)]
#    if (!any(is.na(alleles)) && alleles[1] != alleles[2]) {
#      het_count <- het_count + 1
#    }
#  }
#  het_count
#}

#Calculate Ho across individuals
#ho_individual <- apply(geno, 1, calculate_Ho)

#Combine with population info
#df_ho <- data.frame(
#  Individual = indNames(genind_obj),
#  Population = pop(genind_obj),
#  Ho = ho_individual
#)

#df_ho$Population <- factor(df_ho$Population, levels = c(
#  "Diego_up", "Diego_mid", "Diego_down",
#  "Aguacate_up", "Aguacate_mid", "Aguacate_down",
#  "Paraiso_up", "Paraiso_down"  # no mid for Paraiso
#))
##Summarize average and SD per pop
# Summarize average and standard deviation of Ho per population
#ho_summary <- df_ho %>%
#  group_by(Population) %>%
#  summarise(
#    Mean_Ho = mean(Ho, na.rm = TRUE),
#    SD_Ho = sd(Ho, na.rm = TRUE),
#    N = n()
#  )
