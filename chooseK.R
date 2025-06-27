##Admixture results K values
setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography/")
## Marathrum_Panama
K<-c(0.35580, 0.50033, 0.32892, 0.32851, 0.34281, 0.35786, 0.38121, 0.41667,
     0.43273, 0.47852)
clust<-c(1,2,3,4,5,6,7,8,9,10)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Panama.pdf", width = 6, height = 5)

plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Panama")

dev.off()
#K=3-4

## Marathrum Cocle
K<-c(0.32072, 0.55218, 0.31744, 0.33277, 0.35024, 0.37893, 0.41282,
     0.44657, 0.47381, 0.54541)
clust<-c(1,2,3,4,5,6,7,8,9,10)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Colon.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Colón")
dev.off()
#K=3

## Marathrum Chiriqui
K<-c(0.36086, 0.32108, 0.54336, 0.61848, 0.72954, 0.79105, 0.40839, 0.89348,
     0.62648, 0.92185)
clust<-c(1,2,3,4,5,6,7,8,9,10)
data <- data.frame(clust, K)
data
plot(data)

pdf("cv_error_plot_Marathrum_Chiriqui.pdf", width = 6, height = 5)
plot(data$clust, data$K, type = "p", pch = 19, col = "black",
     xlab = "Number of Clusters (K)", ylab = "Cross-Validation Error",
     main = "Marathrum Chiriquí")
dev.off()
#K=2
