library(vcfR)
library(adegenet)
library(hierfstat)
library(snpR)

setwd("~/Bedoya Dropbox/Bedoya_Research_Group/River_phylogeography")


library(snpR)
library(vcfR)
library(adegenet)
library(hierfstat)

vcf_path  <- "marathrum_panama_filtered_final_bin.recode.vcf"
pops_path <- "pops.txt"

pops <- read.table(pops_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
pops$sample <- trimws(pops$sample)
pops$pop    <- trimws(pops$pop)

############################
## FIS estimation (snpR)
############################
x <- read_vcf(vcf_path)

samps <- trimws(as.character(sample.meta(x)$sampID))
sample.meta(x)$pop <- pops$pop[match(samps, pops$sample)]

stopifnot(!any(is.na(sample.meta(x)$pop)))  # fails if any sample names don't match

x <- calc_fis(x, facets="pop")
fis_pop <- get.snpR.stats(x, facets="pop", stats="fis")
fis_pop


##################
## CI estimation
##################
vcf <- read.vcfR(vcf_path)
gi  <- vcfR2genind(vcf)

pop_vec <- pops$pop[match(indNames(gi), pops$sample)]
stopifnot(!any(is.na(pop_vec)))

pop(gi) <- factor(pop_vec)

dat_hf <- genind2hierfstat(gi)
boot_fis <- boot.ppfis(dat_hf, nboot=9999)
boot_fis

#levels(pop(gi))
#table(pop(gi))

ci <- boot_fis$fis.ci
pop_names <- levels(pop(gi))

out <- data.frame(
  pop = pop_names,
  CI_low = ci[, "ll"],
  CI_high = ci[, "hl"],
  row.names = NULL
)
out

#######################
####FST Estimation ####
#######################

bs <- basic.stats(dat_hf)
fst_global <- bs$overall["Fst"]
fst_global

------
set.seed(1)
loci_cols <- 2:ncol(dat_hf)

boot_vals <- numeric(nboot)
for (b in seq_len(nboot)) {
  resamp <- sample(loci_cols, length(loci_cols), replace = TRUE)
  dat_b <- dat_hf[, c(1, resamp), drop = FALSE]
  boot_vals[b] <- global_fst(dat_b)
}

ci <- quantile(boot_vals, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

cat("Global FST =", fst_global, "\n")
cat("95% CI =", ci[1], "to", ci[2], "\n")
-----
  
fst_pairwise <- pairwise.WCfst(dat_hf) 
fst_pairwise
ci_mat <- boot.ppfst(dat_hf, nboot=9999) 
ci_mat


