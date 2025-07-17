library(ggplot2)
library(patchwork)


#Diego
diego_samples <- 2:52
diego_sites <- c(493,740,894,1003,1084,1149,1201,1246,1284,1316,1344,1370,1393,
            1413,1432,1448,1464,1476,1489,1498,1509,1518,1527,1534,1542,1545,
            1552,1556,1562,1567,1573,1577,1582,1582,1586,1585,1589,1589,1593,
            1583,1587,1581,1584,1570,1572,1557,1560,1475,1477,996,997)	
diego <- data.frame(
  samples = diego_samples, sites = diego_sites
)

#Aguacate
aguacate_samples <- 2:40
aguacate_sites <- c(539,808,975,1092,1180,1249,1305,1351,1391,1425,1456,1483,
                    1508,1526,1546,1563,1579,1591,1604,1612,1623,1627,1637,
                    1645,1653,1655,1662,1658,1664,1655,1661,1654,1658,1643,1647,
                    1545,1548,1097,1099)	
aguacate <- data.frame(
  samples = aguacate_samples, sites = aguacate_sites
)

#Cocle
cocle_samples <- 2:22
cocle_sites <- c(490,734,888,996,1078,1140,1193,1233,1271,1301,1330,1349,1372,
                 1377,1395,1392,1407,1379,1391,1055,1063)
cocle <- data.frame(
  samples = cocle_samples, sites = cocle_sites
)

par(mfrow = c(3, 1)) 

p1<-ggplot(diego, aes(x = samples, y = sites)) +
  geom_point(size = 1) +  # Adjust point size here
  xlab("Samples in the projection") +  # Customize x-axis label
  ylab("Segregating sites") +  # Customize y-axis label
  theme(axis.text = element_text(size = 10)) # Adjust axis label font size
p2<-ggplot(aguacate, aes(x = samples, y = sites)) +
  geom_point(size = 1) +  # Adjust point size here
  xlab("Samples in the projection") +  # Customize x-axis label
  ylab("Segregating sites") +  # Customize y-axis label
  theme(axis.text = element_text(size = 10)) # Adjust axis label font size             
p3<-ggplot(cocle, aes(x = samples, y = sites)) +
  geom_point(size = 1) +  # Adjust point size here
  xlab("Samples in the projection") +  # Customize x-axis label
  ylab("Segregating sites") +  # Customize y-axis label
  theme(axis.text = element_text(size = 10)) # Adjust axis label font size

combined_plot <- p1 / p2 / p3
combined_plot
