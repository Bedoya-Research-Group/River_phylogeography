##Admixture results K values
setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/")
## Marathrum_Panama
K<-c(0.31935, 0.25531, 0.23834, 0.24238, 0.25139, 0.27305)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Panama.pdf", width = 6, height = 5)

plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Panama")

dev.off()
#K=3

## Marathrum Colon
K<-c(0.24466, 0.22261, 0.23011, 0.23899, 0.24901, 0.26789)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Colon.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Colón")
dev.off()
#K=2

## Marathrum Chiriqui
K<-c(0.26318, 0.35673, 0.37889, 0.46440, 0.41611, 0.28402)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Chiriqui.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Chiriquí")
dev.off()
#K=1, 6 (higher CV)

## Marathrum Diego
K<-c(0.20114, 0.22998, 0.28534, 0.33538, 0.38091, 0.45274)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Diego.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Diego River")
dev.off()
#K=1

## Marathrum Aguacate
K<-c(0.26805, 0.28551, 0.35614, 0.44430, 0.47682, 0.58333)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Aguacate.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Aguacate River")
dev.off()
#K=1-2

## Marathrum Cocle Norte
K<-c(0.26820, 0.39696, 0.53535, 0.62703, 0.57799, 0.62882)
clust<-c(1,2,3,4,5,6)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Cocle_Norte.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Cocle Norte River")
dev.off()
#K=1
