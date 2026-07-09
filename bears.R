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
library(grid)
library(ggplot2)
library(glue) # For convenient definition of a large phylogeny.

###################################################################################################
### PRE-COMPUTE f2-STATISTICS #####################################################################
###################################################################################################

# The data are from the 2023 article:
# "Range-wide whole-genome resequencing of the brown bear reveals drivers of intraspecies divergence" by Menno J. de Jong et al.
# Files are stored at https://doi.org/10.5061/dryad.qbzkh18n6.
# The thinned autosomal DNA data are in Brown135_auto.mysnps.thinned.20000.vcf.gz.
# Allocation into 27 populations are in brown135_popfile.allinfo.txt.
# The first step is to turn the vcf into binary PLINK files:
# plink --vcf Brown135_auto.mysnps.thinned.20000.vcf.gz --allow-extra-chr --make-bed --out bears
fam <- read.table("bears.fam", stringsAsFactors = FALSE)
inds <- fam$V2
popinfo <- read.table("Brown135_popfile.allinfo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dict <- setNames(popinfo[[2]], popinfo[[1]])
pops <- dict[inds]
pops[93] <- "NorthScand"; pops[94] <- "NorthScand" # Assume "Kola1" = "Russia_Kola1" and "Kola3" = "Russia_Kola3".
extract_f2("bears", "bears", overwrite = T, poly_only = T, auto_only = F, inds = inds, pops = pops)
blocks <- f2_from_precomp("bears")
dim(blocks) # There are 2667 blocks of size 27x27.
count_snps(blocks) # There are 90928 variants.

###################################################################################################
### THE PHYLOGENY #################################################################################
###################################################################################################

# I'm using Newick notation with leaf labels only, no branch lengths.
# The tree is from Figure 4 panel g of the de Jong et al. article, the "backbone" tree made by TreeMix.
# For convenience I'm using the R package glue to define the tree.
Europe <- "(Ural,((Baltic,(MiddleEast,Europe)),(NorthScand,(MidScand,SouthScand))))"
Asia <- "((CentreRus,CentreRus2),(Kamtchatka,(Yakutia,(Magadan,(Amur,(Sakhalin,Hokkaido))))))"
America <- "(Kodiak,(Aleutian,(Alaska,(HudsonBay,(ABCcoast1,(ABCbc,(ABCa,(ABCcoast2,Westcoast))))))))"
tree <- glue("(Black,(polar,(Himalaya,({Europe},({Asia},{America})))))")
less_Americans <- glue("((((((Kodiak,(Aleutian,(Alaska,(ABCbc,HudsonBay)))),{Asia}),{Europe}),Himalaya),polar),Black)")

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
### RUNNING MALINSKY'S BRANCH STATISTIC ###########################################################
###################################################################################################

# Computing the statistics:
table <- all_fbranches(blocks, tree, type = "fb", save = "bears_Malinsky")

# Testing the default plotter (not used in the article):
plot_fbranch(table, "bears_Malinsky_default", cellsize = 5, inner = 5, plussize = 3, terminal = 25)

# (Supplementary) figure S3 (heat map only, trees made with LaTeX tikzpicture):
max(table, na.rm = TRUE) # 0.7205772
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
ggsave("bears_Malinsky.pdf", height = 10.731, width = 11.19, units = "cm")

###################################################################################################
### USING ENHANCED BRANCH STATISTICS ##############################################################
###################################################################################################

# Computing the statistics:
table <- all_fbranches(blocks, tree, demand_n_at_least_two = TRUE, save = "bears_enhanced")

# Testing the default plotter (not used in the article):
plot_fbranch(table, "bears_enhanced_default", cellsize = 5, inner = 5, plussize = 3, terminal = 25)

# Figure 7 (heat map only, trees made with LaTeX tikzpicture):
max(table, na.rm = TRUE) # 0.3692701
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
ggsave("bears_enhanced.pdf", height = 10.731, width = 11.19, units = "cm")

###################################################################################################
### EXAMINING THE POLAR BEAR COLUMN MORE CLOSELY ##################################################
###################################################################################################

# The enhanced statistic doesn't detect gene flow from polar into ABCbc, while it's generally recognized as a true event. Why?
# Let's dissect the statistic.

fB(blocks, tree, "polar", "ABCbc", details = TRUE)
# The lower medians for individual (D, E)-pairs:
#
# $result
# [1] 0
#
# $oversets
#       [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]       [,8]       [,9]      [,10]
# [1,]    NA 0.03352305 0.03367351 0.02088248 0.02045607 0.01711295 0.01933921 0.02954044 0.04272552 0.02927579
# [2,]    NA         NA 0.02504419 0.00000000 0.00000000 0.00000000 0.00000000 0.02237737 0.06299100 0.03115674
# [3,]    NA         NA         NA 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.10830993 0.00000000
# [4,]    NA         NA         NA         NA 0.00000000 0.00000000 0.00000000 0.10064674 0.21140978 0.10668272
# [5,]    NA         NA         NA         NA         NA 0.00000000 0.00000000 0.11898984 0.23788319 0.10508761
# [6,]    NA         NA         NA         NA         NA         NA 0.60885891 0.53659497 0.59314217 0.25905166
# [7,]    NA         NA         NA         NA         NA         NA         NA 0.45516880 0.58465656 0.23225139
# [8,]    NA         NA         NA         NA         NA         NA         NA         NA 0.70378073 0.00000000
# [9,]    NA         NA         NA         NA         NA         NA         NA         NA         NA 0.00000000
# [10,]   NA         NA         NA         NA         NA         NA         NA         NA         NA         NA
#
# $sets
# [1]  "polar"                                                                 
# [2]  "Black"                                                                 
# [3]  "Himalaya"                                                              
# [4]  "Baltic,Europe,MiddleEast,MidScand,NorthScand,SouthScand,Ural"          
# [5]  "Amur,CentreRus,CentreRus2,Hokkaido,Kamtchatka,Magadan,Sakhalin,Yakutia"
# [6]  "Kodiak"                                                                
# [7]  "Aleutian"                                                              
# [8]  "Alaska"                                                                
# [9]  "HudsonBay"                                                             
# [10] "ABCcoast1"                                                             
#
# $medians
# [1] 0.03352305 0.02504419 0.00000000 0.00000000 0.00000000 0.00000000 0.10064674 0.21140978 0.03115674
# Certain E return enough negative ratios to collapse the statistic into zero.
# Gene flow from Kodiak to ABCbc is a problem, but likely the biggest issue is too dense sampling among the American brown bears,
# many of which are involved in admixture with polar bears and each other.
# Discarding ABCcoast1, ABCa, ABCcoast2 and Westcoast clarifies the picture.
 
fB(blocks, less_Americans, "polar", "ABCbc", details = TRUE)
# The lower medians for individual (D, E)-pairs:
#
# $result
# [1] 0.04780453
#
# $oversets
#      [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]       [,8]
# [1,]   NA 0.04780453 0.05032491 0.04908311 0.04880464 0.04986558 0.04959242 0.05147107
# [2,]   NA         NA 0.05929640 0.05271410 0.05152677 0.05479069 0.05378483 0.05981026
# [3,]   NA         NA         NA 0.00000000 0.00000000 0.04558783 0.04314495 0.06071967
# [4,]   NA         NA         NA         NA 0.08737836 0.00000000 0.05885747 0.08835331
# [5,]   NA         NA         NA         NA         NA 0.00000000 0.00000000 0.08563491
# [6,]   NA         NA         NA         NA         NA         NA 0.00000000 0.00000000
# [7,]   NA         NA         NA         NA         NA         NA         NA 0.25422414
# [8,]   NA         NA         NA         NA         NA         NA         NA         NA
#
# $sets
# [1] "polar"                                                                 
# [2] "Black"                                                                 
# [3] "Himalaya"                                                              
# [4] "Baltic,Europe,MiddleEast,MidScand,NorthScand,SouthScand,Ural"          
# [5] "Amur,CentreRus,CentreRus2,Hokkaido,Kamtchatka,Magadan,Sakhalin,Yakutia"
# [6] "Kodiak"                                                                
# [7] "Aleutian"                                                              
# [8] "Alaska"                                                                
#
# $medians
# [1] 0.04780453 0.05032491 0.04908311 0.04880464 0.04558783 0.04314495 0.06071967

###################################################################################################
### UNIDIRECTIONAL ADMIXTURE PROPORTION ESTIMATES FOR MANY CHOICES OF D AND E #####################
###################################################################################################

# This is for Figure 8.

B_to_P <- numeric(0)
P_to_B <- numeric(0)
Eurasia <- c("Hokkaido", "Sakhalin", "Amur", "Magadan", "Yakutia", "Kamtchatka", "CentreRus2", "CentreRus",
             "SouthScand", "MidScand", "NorthScand", "Europe", "MiddleEast", "Baltic", "Ural", "Himalaya")
America <- c("Kodiak", "Aleutian", "Alaska", "HudsonBay", "ABCcoast1", "ABCa", "ABCcoast2", "Westcoast")
i <- 1
for (e in Eurasia) {
  for (a in America) {
    B_to_P[i] <- qpdstat(blocks, "polar", "Black", a, e)$est / qpdstat(blocks, "ABCbc", "Black", a, e)$est
    P_to_B[i] <- qpdstat(blocks, "ABCbc", a, e, "Black")$est / qpdstat(blocks, "polar", a, e, "Black")$est
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
    legend.box.spacing = unit(0.275, "cm")
  )
ggsave("B_to_P.pdf", height = 3.247, width = 2.19, units = "cm")

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
    legend.box.spacing = unit(0.275, "cm")
  )
ggsave("P_to_B.pdf", height = 3.247, width = 2.19, units = "cm")

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
unidirectional_alpha(blocks, "Black", "polar", "ABCbc", "Kodiak", "Hokkaido")
# We get alpha = -0.05409189, standard error = 0.07857384.

# Unidirectional estimate of gene flow from polar into ABCbc when A = Westcoast and E = Kamtchatka:
unidirectional_alpha(blocks, "Westcoast", "ABCbc", "polar", "Black", "Kamtchatka")
# We get alpha = -0.1194174, standard error = 0.0140731.

# Bidirectional estimate of gene flow from ABCbc into polar when D = Aleutian and E = Europe:
bidirectional_alpha(blocks, "Black", "polar", "ABCbc", "Aleutian", "Europe")
# We get alpha = 0.2427674, standard error = 0.03439345.

# Unidirectional estimate of gene flow from ABCbc into polar when D = Aleutian and E = Europe:
unidirectional_alpha(blocks, "Black", "polar", "ABCbc", "Aleutian", "Europe")
# We get alpha = 0.2711622, standard error = 0.03727932.

# Bidirectional estimate of gene flow from polar into ABCbc when A = Aleutian and E = Europe:
bidirectional_alpha(blocks, "Aleutian", "ABCbc", "polar", "Black", "Europe")
# We get alpha = 0.1047153, standard error = 0.009286217.

# Unidirectional estimate of gene flow from polar into ABCbc when A = Aleutian and E = Europe:
unidirectional_alpha(blocks, "Aleutian", "ABCbc", "polar", "Black", "Europe")
# We get alpha = 0.1382869, standard error = 0.0110793.

# Bidirectional estimate of gene flow from Westcoast into polar when D = Aleutian and E = Europe:
bidirectional_alpha(blocks, "Black", "polar", "Westcoast", "Aleutian", "Europe")
# We get alpha = 0.2600125, standard error = 0.0373217.

# Unidirectional estimate of gene flow from Westcoast into polar when D = Aleutian and E = Europe:
unidirectional_alpha(blocks, "Black", "polar", "Westcoast", "Aleutian", "Europe")
# We get alpha = 0.3037591, standard error = 0.04129335.

# Bidirectional estimate of gene flow from polar into Westcoast when A = Aleutian and E = Europe:
bidirectional_alpha(blocks, "Aleutian", "Westcoast", "polar", "Black", "Europe")
# We get alpha = 0.1440177, standard error = 0.01106957.

# Unidirectional estimate of gene flow from polar into Westcoast when A = Aleutian and E = Europe:
unidirectional_alpha(blocks, "Aleutian", "Westcoast", "polar", "Black", "Europe")
# We get alpha = 0.1946219, standard error = 0.01177941.

# Like Figure 8 already suggested, looks like Westcoast is more involved with polar bear admixture than the ABC bears are.
qpdstat(blocks, "ABCbc", "Europe", "polar", "Black") # f4 = 0.00672, standard error = 0.000305 
qpdstat(blocks, "Westcoast", "Europe", "polar", "Black") # f4 = 0.00346, standard error = 0.000303 
# But the bare f4-statistic might be emphasizing recent events because branch lengths matter.
