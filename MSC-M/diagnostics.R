# Load the coda package
#install.packages("coda")  # If not already installed
library(coda)
library(knitr)
library(dplyr)
setwd("~/Dropbox/MS/MSCM/Diego/")

# 1. Load the MCMC data (assuming the data is in tabular format)
mcmc1 <- read.table("./up_mid/mcmc1.txt", header = TRUE)
mcmc2 <- read.table("./up_mid/mcmc2.txt", header = TRUE)
mcmc3 <- read.table("./up_mid/mcmc3.txt", header = TRUE)
mcmc4 <- read.table("./up_mid/mcmc4.txt", header = TRUE)

# 2. Combine the chains (columns) into mcmc objects, one for each parameter/column
mcmc_chains <- list()

# For each column (each chain for each parameter), convert to mcmc object
for (i in 1:ncol(mcmc1)) {
  chains <- list(as.mcmc(mcmc1[, i]), 
                 as.mcmc(mcmc2[, i]), 
                 as.mcmc(mcmc3[, i]), 
                 as.mcmc(mcmc4[, i]))
  mcmc_chains[[i]] <- as.mcmc.list(chains)}

# 4. Run Gelman-Rubin diagnostic and ESS for each parameter
gelman_results <- lapply(mcmc_chains, function(chain) {
  result <- gelman.diag(chain)
  result$psrf[, 1]  # Extract the R-hat values (first column)
})

ess_results <- lapply(mcmc_chains, effectiveSize)

# 5. Prepare results for printing
results <- data.frame(
  Parameter = colnames(mcmc1),  # Use column names from the MCMC files
  Gelman_Rhat = sapply(gelman_results, function(x) round(x[1], 2)),  # R-hat values (rounded)
  ESS = sapply(ess_results, function(x) round(x[1], 2))  # Effective Sample Size, rounded
)

# 6. Print the results in a pretty table using knitr::kable
kable(results, format = "markdown", caption = "Gelman-Rubin Diagnostic and ESS for MCMC Parameters")


#7. summarize output
combined_chains <- mcmc.list(
  as.mcmc(mcmc1),
  as.mcmc(mcmc2),
  as.mcmc(mcmc3),
  as.mcmc(mcmc4)
)
combined_mcmc <- as.mcmc(do.call(rbind, combined_chains))
# Calculate mean for each parameter
means <- colMeans(combined_mcmc)

# Calculate HPD intervals
hpd <- HPDinterval(combined_mcmc, prob = 0.95)

# Combine into a single data frame
summary_df <- data.frame(
  Parameter = names(means),
  Mean = round(means, 4),
  `HPD Lower` = round(hpd[, 1], 4),
  `HPD Upper` = round(hpd[, 2], 4),
  row.names = NULL
)

kable(summary_df, caption = "Posterior Summary with 95% HPD Intervals")



