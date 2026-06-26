###################################################################################################
### TESTING AN ENHANCED VERSION OF THE BRANCH STATISTIC ###########################################
###################################################################################################

# Testing an enhanced version of the branch statistic for "Notes on f4-ratio estimation", Kalle Leppälä,
# using the same simulated data set the original was bench marked with.
setwd("C:/Users/kalle/Desktop/cichlids") # Replace with your own folder.
library(admixtools)
library(ggplot2)

###################################################################################################
### PRE-COMPUTE f2-STATISTICS #####################################################################
###################################################################################################

# The simulated data with_geneflow.vcf.gz downloaded from https://github.com/millanek/Dsuite?tab=readme-ov-file.
# Converting into binary PLINK files:
# ./plink --vcf with_geneflow.vcf.gz --make-bed --out simulation
fam <- read.table("simulation.fam", stringsAsFactors = FALSE)
inds <- fam$V2
pops = c(rep(c("m", "n", "o", "p", "q", "r", "s", "t", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), each = 2), "U")
extract_f2("simulation", "simulation", overwrite = T, poly_only = T, inds = inds, pops = pops)
blocks <- f2_from_precomp("simulation")
dim(blocks) # There are 20 blocks.
count_snps(blocks) # There are 4342771 variants.

###################################################################################################
### THE PHYLOGENY #################################################################################
###################################################################################################

# The phylogenetic tree is the file simulated_tree_with_geneflow.nwk from https://github.com/millanek/Dsuite?tab=readme-ov-file.
# The leaf labels have been changed to correspond to Figure 3 in the DSUITE article.
# I'm going to be using a simplified Newick notation that doesn't have branch lengths or commas,
# and assumes that the population label is one symbol. Let the outgroup be "U" (not drawn in the DSUITE article Figure 3).
tree <- "(((((m(n(o(pq))))((r(st))a))(b(cd)))((e(f(g(hi))))(j(kl))))U)"

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
        overB[length(overB) + 1] <- qpdstat(blocks, b, a, source, out)$est / qpdstat(blocks, source, a, source, out)$est 
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
              if (qpdstat(blocks, b, a, d, e)$z > threshold) {
                overB[length(overB) + 1] <- qpdstat(blocks, b, a, d, e)$est / qpdstat(blocks, c, a, d, e)$est
              } else {overB[length(overB) + 1] <- 0}
            }
            overACDE[length(overACDE) + 1] <- min(overB)
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

###################################################################################################
### DSUITE'S FIGURE 3 USING MALINSKY'S F-BRANCH (REPRODUCTION) ####################################
###################################################################################################

rownames <- c("mnopqrstabcd", "mnopqrsta", "mnopq", "m", "nopq", "n", "opq", "o", "pq", "p", "q",
              "rsta", "rst", "r", "st", "s", "t", "a", "bcd", "b", "cd", "c", "d",
              "efghijkl", "efghi", "e", "fghi", "f", "ghi", "g", "hi", "h", "i", "jkl", "j", "kl", "k", "l")
colnames <- c("m", "n", "o", "p", "q", "r", "s", "t", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    table[i, j] <- fb(blocks, tree, j, i, "U")
  }
}
max(table, na.rm = TRUE) # 0.1493714
save(table, file = "fish_Malinsky.RData")
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
    legend.box.spacing = unit(0.2, "cm")
  )
ggsave("fish_Malinsky.pdf", height = 7.6, width = 4.5, units = "cm")

###################################################################################################
### DSUITE'S FIGURE 3 USING ENHANCED F-BRANCH #####################################################
###################################################################################################

rownames <- c("mnopqrstabcd", "mnopqrsta", "mnopq", "m", "nopq", "n", "opq", "o", "pq", "p", "q",
              "rsta", "rst", "r", "st", "s", "t", "a", "bcd", "b", "cd", "c", "d",
              "efghijkl", "efghi", "e", "fghi", "f", "ghi", "g", "hi", "h", "i", "jkl", "j", "kl", "k", "l")
colnames <- c("m", "n", "o", "p", "q", "r", "s", "t", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l")
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    print(c(j, i))
    table[i, j] <- fB(blocks, tree, j, i)
  }
}
max(table, na.rm = TRUE) # 0.1422507
save(table, file = "fish_enhanced.RData")
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
    legend.box.spacing = unit(0.2, "cm")
  )
ggsave("fish_enhanced.pdf", height = 7.6, width = 4.5, units = "cm")

###################################################################################################
### FULL FIGURE USING USING MALINSKY'S F-BRANCH ###################################################
###################################################################################################

rownames <- c("mnopqrstabcdefghijklU", "U", "mnopqrstabcdefghijkl",
              "mnopqrstabcd", "mnopqrsta", "mnopq", "m", "nopq", "n", "opq", "o", "pq", "p", "q",
              "rsta", "rst", "r", "st", "s", "t", "a", "bcd", "b", "cd", "c", "d",
              "efghijkl", "efghi", "e", "fghi", "f", "ghi", "g", "hi", "h", "i", "jkl", "j", "kl", "k", "l")
colnames <- rownames
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    table[i, j] <- fb(blocks, tree, j, i, "U")
  }
}
max(table, na.rm = TRUE) # 0.1493714
save(table, file = "fish_Malinsky_full.RData")
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
ggsave("fish_Malinsky_full.pdf", height = 8.301, width = 8.66, units = "cm")

###################################################################################################
### FULL FIGURE USING USING ENHANCED F-BRANCH #####################################################
###################################################################################################

rownames <- c("mnopqrstabcdefghijklU", "U", "mnopqrstabcdefghijkl",
              "mnopqrstabcd", "mnopqrsta", "mnopq", "m", "nopq", "n", "opq", "o", "pq", "p", "q",
              "rsta", "rst", "r", "st", "s", "t", "a", "bcd", "b", "cd", "c", "d",
              "efghijkl", "efghi", "e", "fghi", "f", "ghi", "g", "hi", "h", "i", "jkl", "j", "kl", "k", "l")
colnames <- rownames
table <- matrix(0, length(rownames), length(colnames))
rownames(table) <- rownames
colnames(table) <- colnames
for (i in rownames) {
  for (j in colnames) {
    print(c(j, i))
    table[i, j] <- fB(blocks, tree, j, i)
  }
}
max(table, na.rm = TRUE) # 0.1919358
save(table, file = "fish_enhanced_full.RData")
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
ggsave("fish_enhanced_full.pdf", height = 8.301, width = 8.66, units = "cm")

###################################################################################################
### TESTING HOW ADMIXTOOL'S FIND_GRAPHS() PERFORMS ################################################
###################################################################################################




