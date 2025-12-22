###################################################################################################
### A REAL LIFE BEAR EXAMPLE ######################################################################
###################################################################################################

# This log describes a light weight bear admixture analysis for "Notes on f4-ratio estimation", Kalle Leppälä.
# The data are from an article "Range-wide whole-genome resequencing of the brown bear reveals drivers of
# intraspecies divergence" by Menno J. de Jong et al, 2023.
# The paper compares signals of structure and gene flow inferred from autosomes, sex chromosomes and mitochondria,
# and discusses how those signals are molded my sex-biased migration and age of events.
# My focus is narrower: I only use the autosomes, aiming to demonstrate the utility of the unrooted
# perspective in f4-ratio estimation and the newly developed enhanced f-branch statistic.
# The main point is that the detection of gene flow from brown bears into polar bears becomes possible,
# while artificial restrictions in f4-ratio on the tree marked with an asterisk or Malinsky's branch statistic prevent it.
setwd("C:/Users/kalle/Desktop/bears") # Replace with your own folder.
library(admixtools)
library(ggplot2)

###################################################################################################
### PRE-COMPUTE f2-STATISTICS #####################################################################
###################################################################################################

# The data are from the 2023 article:
# "Range-wide whole-genome resequencing of the brown bear reveals drivers of intraspecies divergence" by Menno J. de Jong et al.
# Files are stored at https://doi.org/10.5061/dryad.qbzkh18n6.
# The thinned autosomal DNA data are in Brown135_auto.mysnps.thinned.20000.vcf.gz.
# Allocation into 27 populations are in brown135_popfile.allinfo.txt.
# But my code works on population names of one symbol only.
# The English language only has 26 letters so I will assign "0" to Black, and let the other populations be coded as:
# "a": MiddleEast
# "b": Himalaya
# "c": Europe
# "d": SouthScand
# "e": MidScand
# "f": NorthScand
# "g": Baltic
# "h": Ural
# "i": CentreRus2
# "j": CentreRus
# "k": Yakutia
# "l": Amur
# "m": Hokkaido
# "n": Sakhalin
# "o": Magadan
# "p": Kamtchatka
# "q": Aleutian
# "r": Kodiak
# "s": Alaska
# "t": ABCa
# "u": ABCbc
# "v": ABCcoast2
# "w": Westcoast
# "x": ABCcoast1
# "y": HudsonBay
# "z": polar
# The first step is to turn the vcf into binary PLINK files:
# plink --vcf Brown135_auto.mysnps.thinned.20000.vcf.gz --allow-extra-chr --make-bed --out bears
fam <- read.table("bears.fam", stringsAsFactors = FALSE)
inds <- fam$V2
popinfo <- read.table("Brown135_popfile.allinfo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dict1 <- setNames(popinfo[[2]], popinfo[[1]])
popletters <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n",
                "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "0")
popnames <- c("MiddleEast", "Himalaya", "Europe", "SouthScand", "MidScand", "NorthScand", "Baltic",
              "Ural", "CentreRus2", "CentreRus", "Yakutia", "Amur", "Hokkaido", "Sakhalin",
              "Magadan", "Kamtchatka", "Aleutian", "Kodiak", "Alaska", "ABCa", "ABCbc",
              "ABCcoast2", "Westcoast", "ABCcoast1", "HudsonBay", "polar", "Black")
dict2 <- setNames(popletters, popnames)
pops <- dict2[dict1[inds]]
pops[93] <- "f"; pops[94] <- "f" # Assume "Kola1" = "Russia_Kola1" and "Kola3" = "Russia_Kola3".
extract_f2("bears", "bears", overwrite = T, poly_only = T, auto_only = F, inds = inds, pops = pops)
blocks <- f2_from_precomp("bears")
dim(blocks) # There are 2667 blocks of size 27x27.
count_snps(blocks) # There are 90928 variants.

###################################################################################################
### THE PHYLOGENY #################################################################################
###################################################################################################

# I'm going to be using a simplified Newick notation that doesn't have branch lengths or commas,
# and assumes that the population label is one symbol.
# The tree is from Figure 4 panel g of the de Jong et al. article, the "backbone" tree made by TreeMix.
tree <- "(((((((((((((vw)t)u)x)y)s)q)r)((((((mn)l)o)k)p)(ij)))((((de)f)((ac)g))h))b)z)0)"
# I will also need another tree where the number of American brown bear populations is reduced.
less_Americans <- "(((((((((uy)s)q)r)((((((mn)l)o)k)p)(ij)))((((de)f)((ac)g))h))b)z)0)"

###################################################################################################
### FUNCTIONS #####################################################################################
###################################################################################################

# Remove all parentheses from a string.
strip <- function(string) {gsub("[()]", "", string)}

# Order the characters of a string alphabetically.
canonize <- function(string) {paste(sort(strsplit(string, "")[[1]]), collapse = "")}

# Return all descendants of the sister branch of a given branch.
A <- function(tree, branch) {
  # For comparison, make sure branch is in a canonized form:
  branch <- canonize(strip(branch))
  # Remove outer parentheses:
  wood <- substr(tree, 2, nchar(tree) - 1)
  # Find the deepest branching point of the tree:
  result <- ""
  for (i in seq(1, nchar(wood) - 1)) {
    if (sum(strsplit(substr(wood, 1, i), "")[[1]] == "(") == sum(strsplit(substr(wood, 1, i), "")[[1]] == ")")) {
      left <- substr(wood, 1, i)
      right <- substr(wood, i + 1, nchar(wood))
      # If branch is encountered, return its sister branch.
      if (canonize(strip(left)) == branch) {result <- canonize(strip(right))}
      if (canonize(strip(right)) == branch) {result <- canonize(strip(left))}
      # Unless the two branches are terminal, recursively call A.
      if (nchar(left) > 1) {result <- paste0(result, A(left, branch))}
      if (nchar(right) > 1) {result <- paste0(result, A(right, branch))}
    }
  }
  return(result)
}

# Malinsky's f-branch statistic f_b.
fb <- function(blocks, tree, source, branch, out) {
  A <- strsplit(A(tree, branch), "")[[1]]
  B <- strsplit(branch, "")[[1]]
  # Malinsky's f-branch is undefined in certain edge cases:
  if (source %in% A || source %in% B || nchar(source) > 1 || out %in% A || out %in% B || nchar(out) > 1 || source == out) {
    result <- NA
  } else {
    # Median over elements of A:
    overA <- numeric(0)
    for (a in A) {
      overB <- numeric(0)
      # Minimum over elements of B:
      for(b in B) {
        overB[length(overB) + 1] <- qpdstat(blocks, b, a, source, out)$est/qpdstat(blocks, source, a, source, out)$est 
      }
      overA[length(overA) + 1] <- max(0, min(overB))
    }
    result <- median(overA)
  }
  return(result)
}

# Break a tree given in simplified Newick notation into a vector of branches.
branches <- function(tree) {
  if (nchar(tree) == 1) {return(c(tree))}
  else {
    # Remove outer parentheses:
    wood <- substr(tree, 2, nchar(tree) - 1)
    # Find the deepest branching point of the tree:
    for (i in seq(1, nchar(wood) - 1)) {
      if (sum(strsplit(substr(wood, 1, i), "")[[1]] == "(") == sum(strsplit(substr(wood, 1, i), "")[[1]] == ")")) {
        # The root is a branch, and recursively find the rest of the branches of each side.
        return(c(canonize(strip(wood)), branches(substr(wood, 1, i)), branches(substr(wood, i + 1, nchar(wood)))))
      }
    }
  }
}

# Check if x is a subset of y.
subset <- function(x, y) {
  all(strsplit(x, "")[[1]] %in% strsplit(y, "")[[1]])
}

# Given a tree and source and target branches, return the list of connecting sets C = S_0, S_1, ..., S_n, A.
sets <- function(tree, source, target) {
  # This code is a bit ad hoc, and it's commented using rooted language like "sisters", by which I actually
  # mean something closer to the sequence {"siblings", "cousins", "second cousins", ...}.
  # This is because Newick notation is rooted. Sometimes however the leaves in the set C are not actually
  # descendants of the source branch s, whereas the target branch t is.
  # The language might be therefore misleading, but that case is also properly dealt with.
  all <- branches(tree)
  # All ancestors of the source branch:
  sourceancestors <- rev(all[sapply(all, subset, x = source)])
  # All ancestors of the target branch:
  targetancestors <- rev(all[sapply(all, subset, x = target)])
  # Remove overlap:
  overlap <- intersect(sourceancestors, targetancestors)
  sourceancestors <- setdiff(sourceancestors, overlap)
  targetancestors <- setdiff(targetancestors, overlap)
  lastancestor <- overlap[which.min(nchar(overlap))]
  # Sisters of ancestors of the source branch (which counts):
  if (length(sourceancestors) == 0) {
    sourcesisters <- character(0)
    targetancestors <- c(targetancestors, source)
  } else if (length(sourceancestors) == 1) {
    sourcesisters <- source
  } else {
    sourceancestors <- lapply(sourceancestors, function(s) strsplit(s, "")[[1]])
    sourcesisters <- mapply(setdiff, sourceancestors[- 1], sourceancestors[- length(sourceancestors)], SIMPLIFY = FALSE)
    sourcesisters <- c(source, sapply(sourcesisters, paste0, collapse = ""))
  }
  # Sisters of ancestors of the target branch (which doesn't count):
  if (length(targetancestors) == 0) {
    return(NA)
  } else if (length(targetancestors) == 1) {
    targetsisters <- character(0)
  } else {
    targetancestors <- lapply(targetancestors, function(s) strsplit(s, "")[[1]])
    targetsisters <- mapply(setdiff, targetancestors[- 1], targetancestors[- length(targetancestors)], SIMPLIFY = FALSE)
    targetsisters <- rev(sapply(targetsisters, paste0, collapse = ""))
  }
  # Other descendants of the common ancestors:
  othersisters <- setdiff(strsplit(all[1], "")[[1]], strsplit(lastancestor, "")[[1]])
  othersisters <- paste0(othersisters, collapse = "")
  result <- c(sourcesisters, othersisters, targetsisters)
  result <- result[nzchar(result)]
  return(result)
}

# A median function that for even cardinality chooses the lower middle value instead of their average.
lower_median <- function(x, na.rm = FALSE) {
  if (na.rm) {x <- x[!is.na(x)]}
  x <- sort(x)
  n <- length(x)
  if (n == 0) {return(NA)}
  if (n %% 2 == 1) {return(x[(n + 1) / 2])}
  else {return(x[n / 2])}
}

# An enhanced version of the f-branch statistic f_B.
fB <- function(blocks,
               tree,
               source,
               target,
               demand_n_at_least_two = FALSE, # Many false positives arise when n = 1.
               threshold = 2.326348, # Threshold for Z-score when testing that f4(C,A;D,E) is positive.
               details = FALSE)  { # Switch to TRUE for a more detailed output.
  # All possible sets:
  sets <- sets(tree, source, target)
  # Sets A, B and C:
  A <- sets[length(sets)]
  B <- target
  C <- sets[1]
  # Sets form which to choose sets D and E:
  sets <- sets[1:length(sets) - 1]
  if (length(sets) < 2 || (length(sets) == 2 && demand_n_at_least_two == TRUE) || any(is.na(sets))) {
    result <- NA
    oversets <- NA
    medians <- NA
  } else {
    oversets <- matrix(NA, nrow = length(sets), ncol = length(sets))
    # Looping over possible sets D and E:
    for (i in seq(1, length(sets) - 1)) {
      D <- sets[i]
      for (j in seq(i + 1, length(sets))) {
        E <- sets[j]
        overACDE <- numeric(0)
        # Median over elements of A, C, D and E:
        for (a in strsplit(A, "")[[1]]) {for (c in strsplit(C, "")[[1]]) {for (d in strsplit(D, "")[[1]]) {for (e in strsplit(E, "")[[1]]) {
          # Make sure the denominator is significantly positive:
          if (qpdstat(blocks, c, a, d, e)$z > threshold) {
            overB <- numeric(0)
            # Minimum over elements of B:
            for (b in strsplit(B, "")[[1]]) {
              overB[length(overB) + 1] <- qpdstat(blocks, b, a, d, e)$est / qpdstat(blocks, c, a, d, e)$est
            }
            overACDE[length(overACDE) + 1] <- max(0, min(overB))
          }
        }}}}
        # If all the denominators failed to be positive, this (lower) median might become NA:
        oversets[i, j] <- lower_median(overACDE)
      }
    }
    # For each set E, take the (lower) medians over the sets D:
    medians <- numeric(0)
    for (j in seq(2, length(sets))) {medians[length(medians) + 1] <- lower_median(oversets[, j], na.rm = TRUE)}
    # Minimize over the sets E:
    # If all the denominators everywhere failed to be positive, the minimum of medians is also considered NA:
    if (all(is.na(medians))) {
      result <- NA
    } else {
      result <- min(medians, na.rm = TRUE)
    }
  }
  if (details == TRUE) {return(list(result = result, oversets = oversets, sets = sets, medians = medians))}
  else {return(result)}
}

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
### RUNNING ENHANCED FBRANCH ON THE SETTING OF FIGURE 4 LOWER LEFT PANEL ##########################
###################################################################################################

rownames <- branches(tree)
colnames <- rownames
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    print(c(j, i))
    table[i, j] <- fB(blocks, tree, j, i, demand_n_at_least_two = TRUE)
  }
}
# I have to change the order of rows and columns so that the "backbone" tree looks tidy.
order <- c("0abcdefghijklmnopqrstuvwxyz", "0", "abcdefghijklmnopqrstuvwxyz", "z", "abcdefghijklmnopqrstuvwxy", "b",
           "acdefghijklmnopqrstuvwxy", "acdefgh", "h", "acdefg", "acg", "g", "ac", "a", "c", "def", "f", "de", "e", "d",
           "ijklmnopqrstuvwxy", "ijklmnop", "ij", "j", "i", "klmnop", "p", "klmno", "k", "lmno", "o", "lmn", "l", "mn", "n", "m",
           "qrstuvwxy", "r", "qstuvwxy", "q", "stuvwxy", "s", "tuvwxy", "y", "tuvwx", "x", "tuvw", "u", "tvw", "t", "vw", "v", "w")
table <- table[order, order]
max(table, na.rm = TRUE) # 0.4254846
save(table, file = "bears_enhanced.RData")
df <- data.frame(
  row = rep(rownames(table), each = ncol(table)),
  col = rep(colnames(table), times = nrow(table)),
  value = as.vector(t(table))
)
df$row <- factor(df$row, levels = rev(rownames(table)))
df$col <- factor(df$col, levels = colnames(table))
df$over <- sapply(df$value > 0.1, isTRUE)
df$value <- pmin(df$value, 0.1, na.rm = FALSE)

plot <- ggplot(df, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "black", size = 0.3) +
  geom_text(data = df[df$over == TRUE, ], aes(label = "+"), size = 2) +
  scale_fill_gradient(
    low = "white", high = "#4CBB17",
    na.value = "#CCCCCC",
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.01),
    name = NULL,
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.3,
      barwidth = unit(0.2, "cm"),
      barheight = unit(5, "cm"),
      label = FALSE
    )
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0.1, 0, 0, 0),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm")
  )
ggsave("bears_enhanced.pdf", height = 10.731, width = 11.09, units = "cm")

###################################################################################################
### RUNNING MALINSKY'S FBRANCH ON THE SETTING OF FIGURE 4 LOWER LEFT PANEL ########################
###################################################################################################

rownames <- branches(tree)
colnames <- rownames
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    print(c(j, i))
    table[i, j] <- fb(blocks, tree, j, i, "0")
  }
}
order <- c("0abcdefghijklmnopqrstuvwxyz", "0", "abcdefghijklmnopqrstuvwxyz", "z", "abcdefghijklmnopqrstuvwxy", "b",
           "acdefghijklmnopqrstuvwxy", "acdefgh", "h", "acdefg", "acg", "g", "ac", "a", "c", "def", "f", "de", "e", "d",
           "ijklmnopqrstuvwxy", "ijklmnop", "ij", "j", "i", "klmnop", "p", "klmno", "k", "lmno", "o", "lmn", "l", "mn", "n", "m",
           "qrstuvwxy", "r", "qstuvwxy", "q", "stuvwxy", "s", "tuvwxy", "y", "tuvwx", "x", "tuvw", "u", "tvw", "t", "vw", "v", "w")
table <- table[order, order]
max(table, na.rm = TRUE) # 0.7205772
save(table, file = "bears_Malinsky.RData")
df <- data.frame(
  row = rep(rownames(table), each = ncol(table)),
  col = rep(colnames(table), times = nrow(table)),
  value = as.vector(t(table))
)
df$row <- factor(df$row, levels = rev(rownames(table)))
df$col <- factor(df$col, levels = colnames(table))
df$over <- sapply(df$value > 0.1, isTRUE)
df$value <- pmin(df$value, 0.1, na.rm = FALSE)

plot <- ggplot(df, aes(x = col, y = row, fill = value)) +
  geom_tile(color = "black", size = 0.3) +
  geom_text(data = df[df$over == TRUE, ], aes(label = "+"), size = 2) +
  scale_fill_gradient(
    low = "white", high = "#4CBB17",
    na.value = "#CCCCCC",
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.01),
    name = NULL,
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.3,
      barwidth = unit(0.2, "cm"),
      barheight = unit(5, "cm"),
      label = FALSE
    )
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0.1, 0, 0, 0),
    legend.position = "right",
    legend.box.spacing = unit(0.1, "cm")
  )
ggsave("bears_Malinsky.pdf", height = 10.731, width = 11.09, units = "cm")

###################################################################################################
### EXAMINING THE POLAR BEAR COLUMN MORE CLOSELY ##################################################
###################################################################################################

# The enhanced statistic doesn't detect gene flow from polar ("z") into ABCbc ("u"), while it's generally recognized as a true event. Why?
# Let's dissect the statistic.

fB(blocks, tree, "z", "u", details = TRUE)
# The lower medians for individual (D, E)-pairs:
#
# $result
# [1] 0
#
# $oversets
# [,1]       [,2]       [,3]        [,4]        [,5]       [,6]        [,7]       [,8]       [,9]       [,10]
# [1,]   NA 0.03352305 0.03367351 0.020882479 0.020456068 0.01711295 0.019339214 0.02954044 0.04272552 0.029275785
# [2,]   NA         NA 0.02504419 0.006918781 0.005390131 0.00000000 0.007119573 0.02237737 0.06299100 0.031156743
# [3,]   NA         NA         NA 0.000000000 0.000000000 0.00000000 0.000000000 0.01761465 0.10830993 0.026897414
# [4,]   NA         NA         NA          NA 0.000000000 0.00000000 0.006337683 0.10064674 0.21140978 0.106682725
# [5,]   NA         NA         NA          NA          NA 0.00000000 0.000000000 0.11898984 0.23788319 0.105087610
# [6,]   NA         NA         NA          NA          NA         NA 0.608858906 0.53659497 0.59314217 0.259051655
# [7,]   NA         NA         NA          NA          NA         NA          NA 0.45516880 0.58465656 0.232251390
# [8,]   NA         NA         NA          NA          NA         NA          NA         NA 0.70378073 0.007236471
# [9,]   NA         NA         NA          NA          NA         NA          NA         NA         NA 0.000000000
# [10,]   NA         NA         NA          NA          NA         NA          NA         NA         NA          NA
#
# $sets
# [1] "z"        "0"        "b"        "acdefgh"  "ijklmnop" "r"        "q"        "s"        "y"        "x"       
#
# $medians
# [1] 0.033523047 0.025044189 0.006918781 0.000000000 0.000000000 0.006337683 0.100646745 0.211409780 0.031156743

# Certain E return enough negative rations to collapse the statistic into zero.
# Gene flow from Kodiak to ABCbc is a problem, but likely the biggest issue is too dense sampling among the American brown bears,
# many of which are involved in admixture with polar bears and each other.
# Discarding ABCcoast1, ABCa, ABCcoast2 and Westcoast clarifies the picture.

fB(blocks, less_Americans, "z", "u", details = TRUE)
# The lower medians for individual (D, E)-pairs:
#
# $result
# [1] 0.04780453
#
# $oversets
# [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]       [,8]
# [1,]   NA 0.04780453 0.05032491 0.04908311 0.04880464 0.04986558 0.04959242 0.05147107
# [2,]   NA         NA 0.05929640 0.05271410 0.05152677 0.05479069 0.05378483 0.05981026
# [3,]   NA         NA         NA 0.02674209 0.02628375 0.04558783 0.04314495 0.06071967
# [4,]   NA         NA         NA         NA 0.09330431 0.06581236 0.05885747 0.08835331
# [5,]   NA         NA         NA         NA         NA 0.05585979 0.04870204 0.08563491
# [6,]   NA         NA         NA         NA         NA         NA         NA 0.15889733
# [7,]   NA         NA         NA         NA         NA         NA         NA 0.25422414
# [8,]   NA         NA         NA         NA         NA         NA         NA         NA
#
# $sets
# [1] "z"        "0"        "b"        "acdefgh"  "ijklmnop" "r"        "q"        "s"       
#
# $medians
# [1] 0.04780453 0.05032491 0.04908311 0.04880464 0.05479069 0.04959242 0.08563491

###################################################################################################
### UNIDIRECTIONAL ADMIXTURE PROPORTION ESTIMATES FOR MANY CHOICES OF D AND E #####################
###################################################################################################

B_to_P <- numeric(0)
P_to_B <- numeric(0)
Eurasia <- c("m", "n", "l", "o", "k", "p", "i", "j", "d", "e", "f", "c", "a", "g", "h", "b")
America <- c("r", "q", "s", "y", "x", "t", "v", "w")
i <- 1
for (e in Eurasia) {
  for (a in America) {
    B_to_P[i] <- qpdstat(blocks, "z", "0", a, e)$est / qpdstat(blocks, "u", "0", a, e)$est
    P_to_B[i] <- qpdstat(blocks, "u", a, e, "0")$est / qpdstat(blocks, "z", a, e, "0")$est
    print(paste(e, a, B_to_P[i], P_to_B[i]))
    i <- i + 1
  }
}
df <- data.frame(B_to_P = B_to_P, P_to_B = P_to_B, Eurasia = rep(Eurasia, each = length(America)), America = rep(America, length(Eurasia)))
df$Eurasia <- factor(df$Eurasia, levels = Eurasia)
df$America <- factor(df$America, levels = America)

plot <- ggplot(df, aes(x = America, y = Eurasia, fill = B_to_P)) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#C80815", mid = "white", high = "#4CBB17",
    na.value = "#CCCCCC",
    limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, by = 0.01),
    name = NULL,
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.3,
      barwidth = unit(0.2, "cm"),
      barheight = unit(2.5, "cm"),
      label = FALSE
    )
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0.1, 0, 0, 0),
    legend.position = "right",
    legend.box.spacing = unit(0.3, "cm")
  )
ggsave("B_to_P.pdf", height = 3.247, width = 2.2, units = "cm")

plot <- ggplot(df, aes(x = "", y = B_to_P)) +
  geom_boxplot(
    width = 0.2,
    size = 0.3,
    colour = "black",
    outlier.colour = "black",
    outlier.size = 0.3
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5), labels = NULL, expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x   = element_line(color = "black", size = 0.3, lineend = "square"),
    axis.ticks.x  = element_line(color = "black", size = 0.3, lineend = "square"),
    plot.margin = margin(0, 4, 0, 4, unit = "pt")
  )
ggsave("B_to_P_box.pdf", height = 0.812, width = 1.98, units = "cm")

plot <- ggplot(df, aes(x = America, y = Eurasia, fill = P_to_B)) +
  geom_tile(color = "black", size = 0.3) +
  scale_fill_gradient2(
    low = "#C80815", mid = "white", high = "#4CBB17",
    na.value = "#CCCCCC",
    limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, by = 0.01),
    name = NULL,
    guide = guide_colorbar(
      ticks = TRUE,
      ticks.colour = NA,
      frame.colour = "black",
      frame.linewidth = 0.3,
      barwidth = unit(0.2, "cm"),
      barheight = unit(2.5, "cm"),
      label = FALSE
    )
  ) +
  scale_x_discrete(position = "top") +
  coord_fixed(ratio = 1) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(0.1, 0, 0, 0),
    legend.position = "right",
    legend.box.spacing = unit(0.3, "cm")
  )
ggsave("P_to_B.pdf", height = 3.247, width = 2.2, units = "cm")

plot <- ggplot(df, aes(x = "", y = P_to_B)) +
  geom_boxplot(
    width = 0.2,
    size = 0.3,
    colour = "black",
    outlier.colour = "black",
    outlier.size = 0.3
  ) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = c(-0.5, 0, 0.5), labels = NULL, expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(expand = expansion(add = 0.2)) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x   = element_line(color = "black", size = 0.3, lineend = "square"),
    axis.ticks.x  = element_line(color = "black", size = 0.3, lineend = "square"),
    plot.margin = margin(0, 4, 0, 4, unit = "pt")
  )
ggsave("P_to_B_box.pdf", height = 0.812, width = 1.98, units = "cm")

###################################################################################################
### SELECTED ESTIMATORS OF ADMIXTURE PROPORTION WITH BLOCK JACKKNIFING ############################
###################################################################################################

# Unidirectional estimate of gene flow from ABCbc into polar when D = Kodiak and E = Hokkaido:
unidirectional_alpha(blocks, "0", "z", "u", "r", "m")
# We get alpha = -0.05409189, standard error = 0.07857384.

# Unidirectional estimate of gene flow from polar into ABCbc when A = Westcoast and E = Kamtchatka:
unidirectional_alpha(blocks, "w", "u", "z", "0", "p")
# We get alpha = -0.1194174, standard error = 0.0140731.

# Bidirectional estimate of gene flow from ABCbc into polar when D = Aleutian and E = Europe:
bidirectional_alpha(blocks, "0", "z", "u", "q", "c")
# We get alpha = 0.2427674, standard error = 0.03439345.

# Unidirectional estimate of gene flow from ABCbc into polar when D = Aleutian and E = Europe:
unidirectional_alpha(blocks, "0", "z", "u", "q", "c")
# We get alpha = 0.2711622, standard error = 0.03727932.

# Bidirectional estimate of gene flow from polar into ABCbc when A = Aleutian and E = Europe:
bidirectional_alpha(blocks, "q", "u", "z", "0", "c")
# We get alpha = 0.1047153, standard error = 0.009286217.

# Unidirectional estimate of gene flow from polar into ABCbc when A = Aleutian and E = Europe:
unidirectional_alpha(blocks, "q", "u", "z", "0", "c")
# We get alpha = 0.1382869, standard error = 0.0110793.

# Bidirectional estimate of gene flow from Westcoast into polar when D = Aleutian and E = Europe:
bidirectional_alpha(blocks, "0", "z", "w", "q", "c")
# We get alpha = 0.2600125, standard error = 0.0373217.

# Unidirectional estimate of gene flow from Westcoast into polar when D = Aleutian and E = Europe:
unidirectional_alpha(blocks, "0", "z", "w", "q", "c")
# We get alpha = 0.3037591, standard error = 0.04129335.

# Bidirectional estimate of gene flow from polar into Westcoast when A = Aleutian and E = Europe:
bidirectional_alpha(blocks, "q", "w", "z", "0", "c")
# We get alpha = 0.1440177, standard error = 0.01106957.

# Unidirectional estimate of gene flow from polar into Westcoast when A = Aleutian and E = Europe:
unidirectional_alpha(blocks, "q", "w", "z", "0", "c")
# We get alpha = 0.1946219, standard error = 0.01177941.

# Like Figure 7 already suggested, looks like Westcoast is more involved with polar bear admixture than the ABC bears are.
qpdstat(blocks, "u", "c", "z", "0") # f4 = 0.00672, standard error = 0.000305 
qpdstat(blocks, "w", "c", "z", "0") # f4 = 0.00346, standard error = 0.000303 
# But the bare f4-statistic is emphasizing recent events because branch lengths matter.
