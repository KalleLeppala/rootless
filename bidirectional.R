###################################################################################################
### BIDIRECTIONAL ADMIXTURE PROPORTION ESTIMATION FROM SIMULATED DATA #############################
###################################################################################################

# Simulating a bidirectional gene flow scenario and testing the bidirectional admixture proportion
# estimation for "Notes on f4-ratio estimation", Kalle Leppälä.
# This log records the admixture proportion analysis of the data simulated in bidirectional.slim and bidirectional.py.
setwd("C:/Users/kalle/Desktop/bidirectional") # Replace with your own folder.
library(admixtools)

###################################################################################################
### PRE-COMPUTE f2-STATISTICS #####################################################################
###################################################################################################

fam <- read.table("merged.fam", stringsAsFactors = FALSE)
inds <- fam$V2
pops = rep(c("p1", "p2", "p3", "p4", "p5"), each = 10)
extract_f2("merged", "merged", overwrite = T, poly_only = T, inds = inds, pops = pops)
blocks <- f2_from_precomp("merged") # This folder is not big so I put it in GitHub.
dim(blocks) # There are 500 blocks of size 5x5.
count_snps(blocks) # There are 7703507 variants.

###################################################################################################
### FUNCTIONS #####################################################################################
###################################################################################################

# Computes the f4-ratio of a unidirectional gene flow, inference done with block jackknife estimation of variance.
# Each call of qpdstat() has an internal block jackknife that I don't use for anything.
# Argument boot = 0 is an attempt to reduce the time spent on those unused variance estimates, but it doesn't seem to do much.
# I also experimented with the delta method, but it sometimes gave negative variance estimates, so I decided that the jackknife is better.
# The estimated gene flow proportion is from C to B, without taking the reverse gene flow from B to C into account.
unidirectional_alpha <- function(blocks, A, B, C, D, E) {
  b <- dim(blocks)[3]
  oneout <- numeric(b)
  for (i in seq(1, b)) {
    temp <- blocks[, , -i]
    oneout[i] <- qpdstat(temp, B, A, D, E, boot = 0)$est / qpdstat(temp, C, A, D, E, boot = 0)$est
  }
  est <- mean(oneout)
  se <- sqrt((b - 1)^2 * var(oneout, na.rm = TRUE) / b)
  z <- est / se
  p <- 2 * pnorm(- abs(z))
  return(list(A = A, B = B, C = C, D = D, E = E, est = est, se = se, z = z, p = p))
}

# Computes the f4-ratio of a bidirectional gene flow, inference done with block jackknife estimation of variance.
# Each call of qpdstat() has an internal block jackknife that I don't use for anything.
# Argument boot = 0 is an attempt to reduce the time spent on those unused variance estimates, but it doesn't seem to do much.
# The estimated gene flow proportion is from C to B, taking the reverse gene flow from B to C into account.
bidirectional_alpha <- function(blocks, A, B, C, D, E) {
  b <- dim(blocks)[3]
  oneout <- numeric(b)
  for (i in seq(1, b)) {
    temp <- blocks[, , -i]
    numerator <- qpdstat(temp, B, A, D, E, boot = 0)$est * (qpdstat(temp, C, D, A, E, boot = 0)$est - qpdstat(temp, B, D, A, E, boot = 0)$est)
    denominator <- qpdstat(temp, B, A, D, E, boot = 0)$est * qpdstat(temp, C, D, A, E, boot = 0)$est - qpdstat(temp, C, A, D, E, boot = 0)$est * qpdstat(temp, B, D, A, E)$est
    oneout[i] <- numerator / denominator
  }
  est <- mean(oneout)
  se <- sqrt((b - 1)^2 * var(oneout, na.rm = TRUE) / b)
  z <- est / se
  p <- 2 * pnorm(- abs(z))
  return(list(A = A, B = B, C = C, D = D, E = E, est = est, se = se, z = z, p = p))
}

###################################################################################################
### COMPUTING THE ESTIMATORS OF ADMIXTURE PROPORTION ##############################################
###################################################################################################

# Estimating alpha naively without considering gene flow of the opposite direction:
unidirectional_alpha(blocks, "p4", "p3", "p2", "p1", "p5")
# We get alpha = 0.5042137, standard error = 0.005970725.

# Estimating beta naively without considering gene flow of the opposite direction:
unidirectional_alpha(blocks, "p1", "p2", "p3", "p4", "p5")
# We get beta = 0.5068581, standard error = 0.005848434.

# Estimating alpha conscious of the bidirectional nature of the gene flow event:
bidirectional_alpha(blocks, "p4", "p3", "p2", "p1", "p5")
# We get alpha = 0.3340101, standard error = 0.004003434.

# Estimating beta conscious of the bidirectional nature of the gene flow event:
bidirectional_alpha(blocks, "p1", "p2", "p3", "p4", "p5")
# We get beta = 0.3375623, standard error = 0.003920055.
