###################################################################################################
### TESTING AN ENHANCED VERSION OF THE BRANCH STATISTIC ###########################################
###################################################################################################

# Testing an enhanced version of the branch statistic for "Notes on f4-ratio estimation", Kalle Leppälä,
# using the same simulated data set the original was bench marked with.
setwd("C:/Users/kalle/Desktop/cichlids") # Replace with your own folder.
library(admixtools)
library(grid)
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
# I'm Newick notation with leaf labels only and without branch lengths.
# Let the outgroup be "U" (not drawn in the DSUITE article Figure 3).
tree <- "(U,((((m,(n,(o,(p,q)))),((r,(s,t)),a)),(b,(c,d))),((e,(f,(g,(h,i)))),(j,(k,l)))))"

###################################################################################################
### FUNCTIONS #####################################################################################
###################################################################################################

# Remove all parentheses from a string. Then interpret it as words separated by commas, then organize those words alphabetically.
canonize <- function(string) {
  return(paste(sort(trimws(strsplit(gsub("[()]", "", string), ",")[[1]])), collapse = ","))
}

# Return all descendants of the sister branch of a given branch.
A <- function(tree, branch) {
  # For comparison, make sure branch is in a canonized form:
  branch <- canonize(branch)
  # Remove outer parentheses:
  wood <- substr(tree, 2, nchar(tree) - 1)
  result <- ""
  # Find the deepest branching point of the tree:
  for (i in gregexpr(",", wood)[[1]]) { # Loop through the commas, which we know do exist.
    if (sum(strsplit(substr(wood, 1, i), "")[[1]] == "(") == sum(strsplit(substr(wood, 1, i), "")[[1]] == ")")) {
      left <- substr(wood, 1, i - 1)
      right <- substr(wood, i + 1, nchar(wood))
      # If branch is encountered, return its sister branch.
      if (canonize(left) == branch) {result <- canonize(right)}
      if (canonize(right) == branch) {result <- canonize(left)}
      # Unless the two branches are terminal, recursively call A. The paste0() is here just to pass the result back from subtasks.
      if (grepl(",", left)) {result <- paste0(result, A(left, branch))}
      if (grepl(",", right)) {result <- paste0(result, A(right, branch))}
    }
  }
  return(result)
}

# Malinsky's f-branch statistic f_b.
fb <- function(blocks, tree, source, branch, out) {
  A <- strsplit(A(tree, branch), ",")[[1]]
  B <- strsplit(branch, ",")[[1]]
  # Malinsky's f-branch is undefined in certain edge cases:
  if (source %in% A || source %in% B || grepl(",", source) || out %in% A || out %in% B || grepl(",", out) || source == out) {
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
  if (!grepl(",", tree)) {return(c(tree))} # No commas means that the string is a single leaf.
  else {
    # Remove outer parentheses:
    wood <- substr(tree, 2, nchar(tree) - 1)
    # Find the deepest branching point of the tree:
    for (i in gregexpr(",", wood)[[1]]) { # Loop through the commas, which we know do exist.
      if (sum(strsplit(substr(wood, 1, i), "")[[1]] == "(") == sum(strsplit(substr(wood, 1, i), "")[[1]] == ")")) {
        # The root is a branch, and recursively find the rest of the branches of each side.
        return(c(canonize(wood), branches(substr(wood, 1, i - 1)), branches(substr(wood, i + 1, nchar(wood)))))
      }
    }
  }
}

# Check if x is a subset of y.
subset <- function(x, y) {
  all(strsplit(x, ",")[[1]] %in% strsplit(y, ",")[[1]])
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
  lastancestor <- overlap[which.min(sapply(overlap, function(x) stringr::str_count(x, ",")))]
  # Sisters of ancestors of the source branch (which counts):
  if (length(sourceancestors) == 0) {
    sourcesisters <- character(0)
    targetancestors <- c(targetancestors, source)
  } else if (length(sourceancestors) == 1) {
    sourcesisters <- source
  } else {
    sourceancestors <- lapply(sourceancestors, function(s) strsplit(s, ",")[[1]])
    sourcesisters <- mapply(setdiff, sourceancestors[- 1], sourceancestors[- length(sourceancestors)], SIMPLIFY = FALSE)
    sourcesisters <- c(source, sapply(sourcesisters, paste0, collapse = ","))
  }
  # Sisters of ancestors of the target branch (which doesn't count):
  if (length(targetancestors) == 0) {
    return(NA)
  } else if (length(targetancestors) == 1) {
    targetsisters <- character(0)
  } else {
    targetancestors <- lapply(targetancestors, function(s) strsplit(s, ",")[[1]])
    targetsisters <- mapply(setdiff, targetancestors[- 1], targetancestors[- length(targetancestors)], SIMPLIFY = FALSE)
    targetsisters <- rev(sapply(targetsisters, paste0, collapse = ","))
  }
  # Other descendants of the common ancestors:
  othersisters <- setdiff(strsplit(all[1], ",")[[1]], strsplit(lastancestor, ",")[[1]])
  othersisters <- paste0(othersisters, collapse = ",")
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
               threshold = 1.644854, # Threshold for Z-score when testing that f4(B,A;D,E) and f(C,A;D,E) are positive.
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
        for (a in strsplit(A, ",")[[1]]) {for (c in strsplit(C, ",")[[1]]) {for (d in strsplit(D, ",")[[1]]) {for (e in strsplit(E, ",")[[1]]) {
          # Make sure the denominator is significantly positive:
          if (qpdstat(blocks, c, a, d, e)$z > threshold) {
            overB <- numeric(0)
            # Minimum over elements of B:
            for (b in strsplit(B, ",")[[1]]) {
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

# A function for calculating all the possible branch statistics into a table.
# The table can be directly fed for the visualizing function plot_fbranch().
all_fbranches <- function(blocks, # The blocks as produced by admixtools. 
                          tree, # The backbone tree in Newick notation without branch lengths where leaves are labeled.
                          type = "fB", # Value "fB" is enhanced branch statistic, "fb" Malinsky's branch statistic.
                          demand_n_at_least_two = FALSE, # When type = "fB" this is carried on to the function fB().
                          details = FALSE, # When type = "fB" this is carried on to the function fB().
                          progress = TRUE,
                          save = "" # Prefix for the RData file the table is saved in, if provided.
) {
  branches <- branches(tree)
  table <- matrix(0, length(branches), length(branches))
  rownames(table) <- branches; colnames(table) <- branches
  if (type == "fB") {
    c <- 0
    for (t in branches) {
      for (s in branches) {
        table[t, s] <- fB(blocks, tree, s, t, demand_n_at_least_two = demand_n_at_least_two, details = details)
        c <- c + 1
        if (progress == TRUE) {print(paste(round(100 * c / length(branches)**2, 1), "%", sep = ""))}
      }
    }
  }
  if (type == "fb") {
    # Finding the outgroup.
    outgroup <- ""
    wood <- substr(tree, 2, nchar(tree) - 1)
    start <- sub(",.*$", "", wood)
    if (grepl("[()]", start) == FALSE) {outgroup <- start}
    end <- sub("^.*,(.*)$", "\\1", wood)
    if (grepl("[()]", end) == FALSE) {outgroup <- end}
    if (outgroup == "") {stop("The tree has no outgroup.")}
    c <- 0
    for (t in branches) {
      for (s in branches) {
        table[t, s] <- fb(blocks, tree, s, t, outgroup)
        c <- c + 1
        if (progress == TRUE) {print(paste(round(100 * c / length(branches)**2, 1), "%", sep = ""))}
      }
    }
  }
  if (save != "") {save(table, file = paste(save, ".RData", sep = ""))} # Saving the result on disk.
  return(table)
}

# This function, using the help from two other functions defined below, plots the table of f-branch statistics.
# It's a heatmap visualizing the value of f-branch statistic from a source branch (column) to a target branch (row).
# There's two trees at the left and at the top pointing the location of the source and target branches within the tree.
# It produces a pdf and a 600 dpi png to the current working directory.
# The plotter requires the R packages ggplot2 and grid installed and loaded.
plot_fbranch <- function(table, # A table of branch statistics with row- and colnames matching the output of the function branches().
                         prefix, # A name (string) for the pdf and png files created, without the file type.
                         lowcolor = "white", # The color encoding no indication of gene flow.
                         highcolor = "#4CBB17", # The color encoding potential gene flow.
                         missingcolor = "gray", # The color when statistic is missing, for example when it doesn't exist.
                         upper = 0.1, # The highest value of those branch statistic that will be plotted.
                         plussize = 2, # The size of the plus-symbol indicating values exceeding "upper", zero means do not draw.
                         cellsize = 8, # The width of the (square) cells and the legend bar in millimeters.
                         gridborders = 0.3, # The width of grid lines in millimeters.
                         margins = 1, # The width of the margins that are padding the different plot components.
                         dottedcolor = "gray", # The color of dotted lines pointing at inner branches of the tree.
                         treelinewidth = "3", # The width of the branches of the tree and the dotted lines in pts (inch / 72).
                         inner = 8, # How long the inner branches are in millimeters.
                         terminal = 8, # How long the terminal branches are in millimeters, make sure this is long enough for labels.
                         fontsize = 10, # Font size.
                         fontface = "plain", # Font style, like "italic" or "bold".
                         fontfamily = "sans", # Font family (can point to different fonts depending on stuff).
                         sourcetarget = c("source", "target"), # What the two trees should be called, empty strings omit the labels.
                         legendheight = 0.618 # The height of the legend bar as a proportion of plot height.
) {
  if (!require("ggplot2", quietly = TRUE)) {stop("Package 'ggplot2' is required but not installed or loaded.")}
  if (!require("grid", quietly = TRUE)) {stop("Package 'grid' is required but not installed or loaded.")}
  # Some computations that will be repeated later but necessary already here for determining the plot size:
  branches <- rownames(table)
  depth <- rep(0, length(branches))
  names(depth) <- branches
  branches_parts <- strsplit(branches, ",")
  for (i in seq(1, length(branches))) {
    branch <- branches[i]
    branch_parts <- strsplit(branch, ",")[[1]]
    is_child <- sapply(branches_parts, function(x) all(x %in% branch_parts))
    children <- branches[is_child]
    if (length(children) > 1) {
      rest <- branches[i:length(branches)][!(branches[i:length(branches)] %in% children)]
      if (any(depth[rest] >= depth[branch])) {depth[children] <- depth[children] + 1}
    }
  }
  maxdepth <- max(depth)
  total_width <- unit(maxdepth * inner + terminal + (nrow(table) + 2) * cellsize + 4 * margins, "mm")
  total_height <- unit(maxdepth * inner + terminal + nrow(table) * cellsize + 4 * margins, "mm")
  # Creating the pdf:
  pdf(
    paste(prefix, ".pdf", sep = ""),
    width = total_width / 25.4, # In inches.
    height = total_height / 25.4 # In inches.
  )
  draw_fbranch_plot(table,
                    lowcolor,
                    highcolor,
                    missingcolor,
                    upper,
                    plussize,
                    cellsize,
                    gridborders,
                    margins,
                    dottedcolor,
                    treelinewidth,
                    inner,
                    terminal,
                    fontsize,
                    fontface,
                    fontfamily,
                    sourcetarget,
                    legendheight
  )
  invisible(dev.off())
  # Creating the png:
  png(
    filename = paste(prefix, ".png", sep = ""),
    type = "cairo",
    res = 600,
    width = 600 * total_width / 25.4, # In pixels given 600 dpi.
    height = 600 * total_height / 25.4 # In pixels given 600 dpi.
  )
  draw_fbranch_plot(table,
                    lowcolor,
                    highcolor,
                    missingcolor,
                    upper,
                    plussize,
                    cellsize,
                    gridborders,
                    margins,
                    dottedcolor,
                    treelinewidth,
                    inner,
                    terminal,
                    fontsize,
                    fontface,
                    fontfamily,
                    sourcetarget,
                    legendheight
  )
  invisible(dev.off())
}

# This helper function actually draws the plot.
draw_fbranch_plot <- function(table,
                              lowcolor,
                              highcolor,
                              missingcolor,
                              upper,
                              plussize,
                              cellsize,
                              gridborders,
                              margins,
                              dottedcolor,
                              treelinewidth,
                              inner,
                              terminal,
                              fontsize,
                              fontface,
                              fontfamily,
                              sourcetarget,
                              legendheight
) {
  branches <- rownames(table)
  # Start by working out the tree structure because it affects the dimensions of plot components:
  depth <- rep(0, length(branches)); width <- rep(0, length(branches))
  names(depth) <- branches; names(width) <- branches
  branches_parts <- strsplit(branches, ",")
  for (i in seq(1, length(branches))) {
    branch <- branches[i]
    branch_parts <- strsplit(branch, ",")[[1]]
    is_child <- sapply(branches_parts, function(x) all(x %in% branch_parts))
    children <- branches[is_child]
    if (length(children) > 1) {
      width[branch] <- max(which(is_child))
      rest <- branches[i:length(branches)][!(branches[i:length(branches)] %in% children)]
      if (any(depth[rest] >= depth[branch])) {depth[children] <- depth[children] + 1}
    }
  }
  width[duplicated(width)] <- 0
  maxdepth <- max(depth)
  # Assign the dimensions then:
  common_dimension <- unit(nrow(table) * cellsize + 2 * margins, "mm")
  tree_dimension <- unit(maxdepth * inner + terminal + 2 * margins, "mm")
  legend_dimension <- unit(2 * cellsize, "mm")
  # The full layout is 2 × 3; heatmap at (2, 2), trees at (2, 1) and (1, 2), and legend at (2, 3):
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 2,
    ncol = 3,
    widths = unit.c(tree_dimension, common_dimension, legend_dimension),
    heights = unit.c(tree_dimension, common_dimension)
  )))
  # Top tree:
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  draw_tree(table,
            cellsize,
            margins,
            dottedcolor,
            treelinewidth,
            inner, terminal,
            fontsize,
            fontface,
            fontfamily,
            sourcetarget,
            depth,
            width,
            TRUE
  )
  popViewport()
  # Left tree:
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  draw_tree(table,
            cellsize,
            margins,
            dottedcolor,
            treelinewidth,
            inner, terminal,
            fontsize,
            fontface,
            fontfamily,
            sourcetarget,
            depth,
            width,
            FALSE
  )
  popViewport()
  # Heatmap:
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  df <- data.frame(
    row = rep(rownames(table), each = ncol(table)),
    col = rep(colnames(table), times = nrow(table)),
    value = as.vector(t(table))
  )
  df$row <- factor(df$row, levels = rev(rownames(table)))
  df$col <- factor(df$col, levels = colnames(table))
  df$over <- sapply(df$value > upper, isTRUE) # Where to draw the plus symbol.
  if (plussize == 0) {df$over <- FALSE} # Size zero would actually draw something while we wish to draw nothing.
  df$value <- pmin(df$value, upper, na.rm = FALSE) # Cap the branch statistics at the upper value for visibility.
  edges <- seq(0.5, nrow(table) + 0.5, by = 1)
  heat <- ggplot(df, aes(x = col, y = row, fill = value)) +
    geom_tile(color = "black", size = gridborders, linejoin = "round", lineend = "round") +
    geom_text(data = df[df$over == TRUE, ], aes(label = "+"), size = plussize) +
    scale_fill_gradient(
      low = lowcolor, high = highcolor,
      na.value = missingcolor,
      guide = "none"
    ) +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(plot.margin = margin(margins, margins, margins, margins))
  grid.draw(ggplotGrob(heat))
  popViewport()
  # Legend:
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  position <- (1.5 * cellsize - margins) / (2 * cellsize)
  grid.rect(x = position, y = 0.5, width = 0.5, height = legendheight,
            gp = gpar(fill = linearGradient(colours = c(lowcolor, highcolor)),
                      col = "black", lwd = 2.834645669 * gridborders))
  upper_char <- as.character(upper)
  if (grepl("\\.", upper_char)) {decimals <- nchar(sub(".*\\.", "", upper_char))}
  else {decimals <- 0}
  legend_pad <- margins / convertUnit(common_dimension, "mm", valueOnly = TRUE)
  grid.text(formatC(0, format = "f", digits = decimals), x = position,
            y = (1 - legendheight) / 2 - legend_pad, just = "top",
            gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
  grid.text(upper, x = position, 
            y = (1 + legendheight) / 2 + legend_pad, just = "bottom",
            gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
  popViewport()
  popViewport()
  invisible(NULL)
}

# This helper function draws the two trees in the plot.
draw_tree <- function(table,
                      cellsize,
                      margins,
                      dottedcolor,
                      treelinewidth,
                      inner,
                      terminal,
                      fontsize,
                      fontface,
                      fontfamily,
                      sourcetarget,
                      depth,
                      width,
                      flip) {
  branches <- rownames(table)
  maxdepth <- max(depth)
  pushViewport(viewport(clip = "off")) # Plot elements like text seeping onto the padding is allowed.
  if (flip) {
    vp_inner <- viewport(
      x = 0.5, # Center the tree canvas inside the padding.
      y = 0.5, # Center the tree canvas inside the padding.
      width = unit(nrow(table) * cellsize, "mm"),
      height  = unit(maxdepth * inner + terminal, "mm"),
      xscale = c(0.5, length(branches) + 0.5), # Labels and branch tips are mid-cells relative to the heatmap.
      yscale = c(maxdepth * inner + terminal, 0),
      clip = "off"
    )
  } else {
    vp_inner <- viewport(
      x = 0.5, # Center the tree canvas inside the padding.
      y = 0.5, # Center the tree canvas inside the padding.
      width  = unit(maxdepth * inner + terminal, "mm"),
      height = unit(nrow(table) * cellsize, "mm"),
      xscale = c(0, maxdepth * inner + terminal),
      yscale = c(length(branches) + 0.5, 0.5), # Labels and branch tips are mid-cells relative to the heatmap.
      clip = "off"
    )
  }
  pushViewport(vp_inner)
  if (flip) {
    pad <- convertHeight(unit(margins, "mm"), "native", valueOnly = TRUE)
  } else {
    pad <- convertWidth(unit(margins, "mm"), "native", valueOnly = TRUE)
  }
  for (i in seq(1, length(branches))) {
    branch <- branches[i]
    if (grepl(",", branch)) { # Inner branches, drawn first so dotted lines never cover the tree.
      if (flip) {
        grid.segments(y0 = maxdepth * inner + terminal, x0 = i, y1 = depth[i] * inner, x1 = i,
                      default.units = "native", gp = gpar(col = dottedcolor, lwd = treelinewidth, lty = "dotted"))
      } else {
        grid.segments(x0 = maxdepth * inner + terminal, y0 = i, x1 = depth[i] * inner, y1 = i,
                      default.units = "native", gp = gpar(col = dottedcolor, lwd = treelinewidth, lty = "dotted"))
      }
    } else { # Terminal branches, shortened to make room for the labels.
      if (flip) {
        label_width <- fontsize * convertHeight(stringWidth(branch), "native", valueOnly = TRUE) / 12
        grid.segments(y0 = depth[i] * inner, x0 = i, y1 = maxdepth * inner + terminal + label_width + pad, x1 = i,
                      default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
        grid.text(branch, y = unit(maxdepth * inner + terminal, "native"), x = unit(i, "native"), just = "left",
                  rot = 90, gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
      } else {
        label_width <- fontsize * convertWidth(stringWidth(branch), "native", valueOnly = TRUE) / 12
        grid.segments(x0 = depth[i] * inner, y0 = i, x1 = maxdepth * inner + terminal - label_width - pad, y1 = i,
                      default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
        grid.text(branch, x = unit(maxdepth * inner + terminal, "native"), y = unit(i, "native"), just = "right",
                  gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
      }
    }
  }
  # The root:
  if (flip) {
    grid.segments(y0 = 0, x0 = 0, y1 = 0, x1 = 1,
                  default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
  } else {
    grid.segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1,
                  default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
  }
  for (i in seq(1, length(branches))) {
    if (width[i] > 0) { # Orthogonal branches connecting others.
      if (flip) {
        grid.segments(y0 = depth[i] * inner, x0 = i, y1 = depth[i] * inner, x1 = width[i],
                      default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
      } else {
        grid.segments(x0 = depth[i] * inner, y0 = i, x1 = depth[i] * inner, y1 = width[i],
                      default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
      }
      if (i > 1) {
        if (flip) { # The tips for the orthogonal branches, except the first one.
          grid.segments(y0 = (depth[i] - 1) * inner, x0 = i, y1 = depth[i] * inner, x1 = i,
                        default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))          
        } else {
          grid.segments(x0 = (depth[i] - 1) * inner, y0 = i, x1 = depth[i] * inner, y1 = i,
                        default.units = "native", gp = gpar(col = "black", lwd = treelinewidth))
        }
      }
    }
  }
  # Labeling the tree.
  if (flip) {
    grid.text(sourcetarget[1], x = unit(0, "native"), y = unit(0 - 2 * pad, "native"), just = "right",
              rot = 90, gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
  } else {
    grid.text(sourcetarget[2], x = unit(0 + 2 * pad, "native"), y = unit(0, "native"), just = "left",
              gp = gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily))
  }
  popViewport()
  popViewport()
  invisible(NULL)
}

###################################################################################################
### USING MALISNKY's BRANCH STATISTICS ############################################################
###################################################################################################

# Computing the statistics:
table <- all_fbranches(blocks, tree, type = "fb", save = "fish_Malinsky")

# Testing the default plotter (not used in the article):
plot_fbranch(table, "fish_Malinsky_default", cellsize = 5, inner = 5, plussize = 3)

# Full plot; (supplementary) Figure S1 (heat map only, trees made with LaTeX tikzpicture):
max(table, na.rm = TRUE) # 0.1493714
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
ggsave("fish_Malinsky_full.pdf", height = 8.301, width = 8.76, units = "cm")

# Restricted plot; Figure 6 left panel (reproduction of DSUITE paper's Figure 3, heat map only, trees made with LaTeX tikzpicture):
rows <- rownames(table)[4:length(rownames(table))] # No root, no outgroup, no everything-but-the-outgroup.
columns <- colnames(table)[!grepl(",", colnames(table))]; columns <- columns[2:length(columns)] # Only terminal branches, no outgroup.
restricted_table <- table[rows, columns]
max(restricted_table, na.rm = TRUE) # 0.1493714
df <- data.frame(
  row = rep(rownames(restricted_table), each = ncol(restricted_table)),
  col = rep(colnames(restricted_table), times = nrow(restricted_table)),
  value = as.vector(t(restricted_table))
)
df$row <- factor(df$row, levels = rev(rownames(restricted_table)))
df$col <- factor(df$col, levels = colnames(restricted_table))
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
### USING ENHANCED BRANCH STATISTICS ##############################################################
###################################################################################################

# Computing the statistics:
table <- all_fbranches(blocks, tree, save = "fish_enhanced")

# Testing the default plotter (not used in the article):
plot_fbranch(table, "fish_enhanced_default", cellsize = 5, inner = 5, plussize = 3)

# Full plot; (supplementary) Figure S2 (heat map only, trees made with LaTeX tikzpicture):
max(table, na.rm = TRUE) # 0.1919358
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
ggsave("fish_enhanced_full.pdf", height = 8.301, width = 8.76, units = "cm")

# Restricted plot; Figure 6 right panel (heat map only, trees made with LaTeX tikzpicture):
rows <- rownames(table)[4:length(rownames(table))] # No root, no outgroup, no everything-but-the-outgroup.
columns <- colnames(table)[!grepl(",", colnames(table))]; columns <- columns[2:length(columns)] # Only terminal branches, no outgroup.
restricted_table <- table[rows, columns]
max(restricted_table, na.rm = TRUE) # 0.1422507
df <- data.frame(
  row = rep(rownames(restricted_table), each = ncol(restricted_table)),
  col = rep(colnames(restricted_table), times = nrow(restricted_table)),
  value = as.vector(t(restricted_table))
)
df$row <- factor(df$row, levels = rev(rownames(restricted_table)))
df$col <- factor(df$col, levels = colnames(restricted_table))
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
### TESTING HOW ADMIXTOOL'S FIND_GRAPHS() PERFORMS ################################################
###################################################################################################

# Let's see how well find_graphs() algorithm of admixtools does with the cichlid data set:
set.seed(1)
graphs <- find_graphs(blocks, outpop = "U", numadmix = 5, numgraphs = 25, stop_gen = 10000, stop_gen2 = 60, plusminus_generations = 100)
save(graphs, file = "findgraphs.RData")
png("comparison.png", type = "cairo", width = 4800, height = 1500, res = 600)
plot_graph(graphs$edges[[order(graphs$score)[1]]], title = paste("penalty =", round(sort(graphs$score)[1], 3)))
dev.off()
# The admixture graph is pretty far from the true graph, not displaying same trees either.
