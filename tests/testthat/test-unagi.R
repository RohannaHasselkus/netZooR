# unit-tests for UNAGI
context("test UNAGI functions")
test_that("[UNAGI] SignificantBreadthFirstSearchU() function yields expected results", {
  
  # Construct a starting network, which will be modified.
  # Here, we expect genes A and B to be connected after 1 hop via TF2, genes A and D
  # to be connected after 1 hop via TF3, and genes A and C to be connected after 2
  # hops via TF4.
  startingNetwork <- data.frame(source = c(rep(c("gene1", "gene2", "gene3", "gene4"), 4)),
                                target = c(rep("geneA", 4), rep("geneB", 4), rep("geneC", 4), rep("geneD", 4)),
                           score = c(-3, 3, 5, -5, 0, 4, 0.0005, -0.5, 0, 0.0005, -1, 4, 0.0005, -2, 5, 3))
  addedNoise1 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                        rep("geneB", 100),
                                                        rep("geneC", 100),
                                                        rep("geneD", 100)),
                            score = stats::rnorm(400) / 1000)
  addedNoise2 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  addedNoise3 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  fullNetwork1 <- rbind(startingNetwork, addedNoise1)
  rownames(fullNetwork1) <- paste(fullNetwork1$source, fullNetwork1$target, sep = "__")
  
  # Create a set of networks close to the original.
  fullNetworks <- fullNetwork1
  fullNetworks[,4] <- fullNetwork1$score - 0.00005
  fullNetworks[,5] <- fullNetwork1$score + 0.00005
  
  # Use all combined scores as the null distribution.
  null <- stats::rnorm(n = nrow(fullNetworks) * 3) / 100
  
  # Test that errors are thrown when appropriate.
  expect_error(SignificantBreadthFirstSearchU(networks = subnetwork,
                                             startingNodes = c("blob", "fish"),
                                             nodesToExclude = c()),
               "ERROR: Starting nodes do not overlap with network nodes")
  expect_error(SignificantBreadthFirstSearchU(networks = subnetwork, 
                                             startingNodes = c("gene1", "gene2"),
                                             nodesToExclude = c()),
               "ERROR: Starting nodes do not overlap with network nodes")
  expect_error(SignificantBreadthFirstSearchU(networks = subnetwork, 
                                             startingNodes = c("geneA", "geneB", "geneC"),
                                             nodesToExclude = c("blob", "fish")),
               "ERROR: List of nodes to exclude does not overlap with network nodes")
  expect_error(SignificantBreadthFirstSearchU(networks = subnetwork, 
                                             startingNodes = c("geneA", "geneB", "geneC"),
                                             nodesToExclude = c("geneA", "geneB")),
               "ERROR: Starting nodes cannot overlap with nodes to exclude")
  
  # Ensure that, when starting from genes A, B, and C, we obtain the correct values.
  expect_true(length(setdiff(c("gene2__geneA", "gene2__geneB", "gene3__geneA", "gene4__geneC"),
                      rownames(SignificantBreadthFirstSearchU(networks = subnetwork,
                                                  startingNodes = c("geneA", "geneB", "geneC"),
                                                  nodesToExclude = c())))) == 0)
  
  # Ensure that, when starting from genes gene2, gene3, and gene4 and removing genes A, B, and C, we obtain the correct values.
  expect_true(length(setdiff(c("gene3__geneD", "gene4__geneD"),
                             rownames(SignificantBreadthFirstSearchU(networks = subnetwork,
                                                                    startingNodes = c("gene2", "gene3", "gene4"),
                                                                    nodesToExclude = c("geneA", "geneB", "geneC"))))) == 0)
  
  # Unipartite structure.
  startingNetwork <- data.frame(source = c(rep(c("gene1", "gene2", "gene3", "gene4"), 4)),
                                target = c(rep("geneA", 4), rep("geneB", 4), rep("geneC", 4), rep("geneD", 4)),
                                score = c(-3, 3, 5, -5, 0, 4, 0.0005, -0.5, 0, 0.0005, -1, 4, 0.0005, -2, 5, 3))
  addedNoise1 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                               rep("geneB", 100),
                                                               rep("geneC", 100),
                                                               rep("geneD", 100)),
                            score = stats::rnorm(400) / 1000)
  
  geneANetHop1 <- data.frame(source = c("gene2", "gene3","geneB"), target = c("geneA", "geneA", "geneA"))
  geneBNetHop1 <- data.frame(source = "gene2", target = "geneB")
  geneCNetHop1 <- data.frame(source = c("gene4","geneD"), target = c("geneC", "geneC"))
  geneDNetHop1 <- data.frame(source = c("gene3", "gene4"), target = c("geneD", "geneD"))
  gene2NetHop1 <- data.frame(source = c("geneA", "geneB", "gene3"), target=c("gene2", "gene2", "gene2"))
  geneANetHop2 <- rbind(geneBNetHop1, geneDNetHop1[1,])
  geneBNetHop2 <- geneANetHop1[1,]
  geneCNetHop2 <- geneDNetHop1[2,]
  geneDNetHop2 <- rbind(geneANetHop1[2,], geneCNetHop1)
  geneANetHop3 <- geneDNetHop1[2,]
  geneBNetHop3 <- geneANetHop1[2,]
  geneCNetHop3 <- geneDNetHop1[1,]
  geneDNetHop3 <- geneANetHop1[1,]
  
  subnetworks <- list(geneA = list(geneANetHop1, geneANetHop2, geneANetHop3),
                      geneB = list(geneBNetHop1, geneBNetHop2, geneBNetHop3),
                      geneC = list(geneCNetHop1, geneCNetHop2, geneCNetHop3),
                      geneD = list(geneDNetHop1, geneDNetHop2, geneDNetHop3))
})
#toy network -- Troubleshooting tests 2-5
subnetwork <- startingNetwork[, c("source", "target")]
rownames(subnetwork) <- paste(subnetwork$source, subnetwork$target, sep = "__")
test_that("[UNAGI] FindConnectionsForAllHopCountsU() function yields expected results", {
  
  # Set up the subnetwork from the previous example.
  geneANetHop1 <- data.frame(source = c("geneA", "geneA"), target = c("gene2", "gene3"))
  geneBNetHop1 <- data.frame(source = "geneB", target = "gene2")
  geneCNetHop1 <- data.frame(source = "geneC", target = "gene4")
  geneDNetHop1 <- data.frame(source = c("geneD", "geneD"), target = c("gene3", "gene4"))
  geneANetHop2 <- data.frame(source = c("gene2", "gene3"), target = c("geneB", "geneD"))
  geneBNetHop2 <- data.frame(source = "gene2", target = "geneA")
  geneCNetHop2 <- data.frame(source = "gene4", target = "geneD")
  geneDNetHop2 <- data.frame(source = c("gene3", "gene4"), target = c("geneA", "geneC"))
  geneANetHop3 <- data.frame(source = "geneD", target = "gene4")
  geneBNetHop3 <- data.frame(source = "geneA", target = "gene3")
  geneCNetHop3 <- data.frame(source = "geneD", target = "gene3")
  geneDNetHop3 <- data.frame(source = "geneA", target = "gene2")
  
  subnetworks <- list(geneA = list(geneANetHop1, geneANetHop2, geneANetHop3),
                      geneB = list(geneBNetHop1, geneBNetHop2, geneBNetHop3),
                      geneC = list(geneCNetHop1, geneCNetHop2, geneCNetHop3),
                      geneD = list(geneDNetHop1, geneDNetHop2, geneDNetHop3))
  
  # Obtain the overlaps for 1, 2, and 3 hops.
  result <- FindConnectionsForAllHopCountsU(subnetworks)
  
  if (is.data.frame(result) || is.matrix(result)) {
    expect_setequal(rownames(result),
                    c("gene2__geneA", "gene2__geneB", "gene3__geneA", "gene3__geneD",
                      "gene4__geneC", "gene4__geneD", "geneA__gene3", "geneD__gene3",
                      "geneA__gene2", "geneB__gene2"))
    
  } else {
    stop("Error: FindConnectionsForAllHopCountsU did not return a data frame or matrix.")
  }
  
  # Obtain the overlaps for only the first two hops, for only A and C.
  subnetworks <- list(geneA = list(geneANetHop1, geneANetHop2),
                      geneC = list(geneCNetHop1, geneCNetHop2))
  result <- FindConnectionsForAllHopCountsU(subnetworks)
  
  if (is.data.frame(result) || is.matrix(result)) {
    expect_equal(rownames(result), 
                 c("geneA__gene3", "gene3__geneD", "geneC__gene4", "gene4__geneD"))
  } else {
    stop("Error: FindConnectionsForAllHopCountsU did not return a data frame or matrix.")
  }
})
test_that("[UNAGI] FindSignificantEdgesforhopU() function yields expected results",{
  
  # Construct a starting network, which will be modified.
  # Here, we expect genes A and B to be connected after 1 hop via gene2, genes A and D
  # to be connected after 1 hop via TF3, and genes A and C to be connected after 2
  # hops via TF4.
  startingNetwork <- data.frame(source = c(rep(c("gene1", "gene2", "gene3", "gene4"), 4)),
                                target = c(rep("geneA", 4), rep("geneB", 4), rep("geneC", 4), rep("geneD", 4)),
                                score = c(-3, 3, 5, -5, 0, 4, 0.0005, -0.5, 0, 0.0005, -1, 4, 0.0005, -2, 5, 3))
  addedNoise1 <- data.frame(source = rep(1:100, 4), gene = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  addedNoise2 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  addedNoise3 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  fullNetwork1 <- rbind(startingNetwork, addedNoise1)
  rownames(fullNetwork1) <- paste(fullNetwork1$gene, fullNetwork1$gene, sep = "__")
  
  # Create a set of networks close to the original.
  fullNetworks <- fullNetwork1
  fullNetworks[,4] <- fullNetwork1$score - 0.00005
  fullNetworks[,5] <- fullNetwork1$score + 0.00005
  
  #  Null distribution is a narrower version of the normal distribution.
  null <- rnorm(n = nrow(fullNetworks) * 3) / 100
  
  # Set up the subnetwork from the previous example.
  geneANetHop1 <- data.frame(source = c("gene2", "gene3"), target = c("geneA", "geneA"))
  rownames(geneANetHop1) <- paste(geneANetHop1$gene, geneANetHop1$gene, sep = "__")
  geneBNetHop1 <- data.frame(source = "gene2", target = "geneB")
  rownames(geneBNetHop1) <- paste(geneBNetHop1$gene, geneBNetHop1$gene, sep = "__")
  geneCNetHop1 <- data.frame(source = "gene4", target = "geneC")
  rownames(geneCNetHop1) <- paste(geneCNetHop1$gene, geneCNetHop1$gene, sep = "__")
  geneDNetHop1 <- data.frame(source = c("gene3", "gene4"), target = c("geneD", "geneD"))
  rownames(geneDNetHop1) <- paste(geneDNetHop1$gene, geneDNetHop1$gene, sep = "__")
  geneANetHop2 <- rbind(geneBNetHop1, geneDNetHop1[1,])
  rownames(geneANetHop2) <- paste(geneANetHop2$gene, geneANetHop2$gene, sep = "__")
  geneBNetHop2 <- geneANetHop1[1,]
  rownames(geneBNetHop2) <- paste(geneBNetHop2$gene, geneBNetHop2$gene, sep = "__")
  geneCNetHop2 <- geneDNetHop1[2,]
  rownames(geneCNetHop2) <- paste(geneCNetHop2$gene, geneCNetHop2$gene, sep = "__")
  geneDNetHop2 <- rbind(geneANetHop1[2,], geneCNetHop1)
  rownames(geneDNetHop2) <- paste(geneDNetHop2$gene, geneDNetHop2$gene, sep = "__")
  geneANetHop3 <- geneDNetHop1[2,]
  rownames(geneANetHop3) <- paste(geneANetHop3$gene, geneANetHop3$gene, sep = "__")
  geneBNetHop3 <- geneANetHop1[2,]
  rownames(geneBNetHop3) <- paste(geneBNetHop3$gene, geneBNetHop3$gene, sep = "__")
  geneCNetHop3 <- geneDNetHop1[1,]
  rownames(geneCNetHop3) <- paste(geneCNetHop3$gene, geneCNetHop3$gene, sep = "__")
  geneDNetHop3 <- geneANetHop1[1,]
  rownames(geneDNetHop3) <- paste(geneDNetHop3$gene, geneDNetHop3$gene, sep = "__")
  subnetworksFull <- list(geneA = list(geneANetHop1, geneANetHop2, geneANetHop3),
                      geneB = list(geneBNetHop1, geneBNetHop2, geneBNetHop3),
                      geneC = list(geneCNetHop1, geneCNetHop2, geneCNetHop3),
                      geneD = list(geneDNetHop1, geneDNetHop2, geneDNetHop3))
  
  # Test method with all genes as input.
  sigEdges <- FindSignificantEdgesforhopU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                                         combinedNetwork = subnetwork,
                                         hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[1]]), rownames(sigEdges$geneB[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[2]]), rownames(sigEdges$geneB[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[3]]), rownames(sigEdges$geneB[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneD[[1]]), rownames(sigEdges$geneD[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneD[[2]]), rownames(sigEdges$geneD[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneD[[3]]), rownames(sigEdges$geneD[[3]]))) == 0)
  
  # Test method with only genes A, B, C.
  sigEdges <- FindSignificantEdgesforhopU(geneSet = c("geneA", "geneB", "geneC"),
                                          combinedNetwork = subnetwork, hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[1]]), rownames(sigEdges$geneB[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[2]]), rownames(sigEdges$geneB[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneB[[3]]), rownames(sigEdges$geneB[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
  
  # Test method with only genes A and C.
  sigEdges <- FindSignificantEdgesforhopU(geneSet = c("geneA", "geneC"),
                                          combinedNetwork = subnetwork,
                                          hopConstraint = 3)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[1]]), rownames(sigEdges$geneA[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[2]]), rownames(sigEdges$geneA[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneA[[3]]), rownames(sigEdges$geneA[[3]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[1]]), rownames(sigEdges$geneC[[1]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[2]]), rownames(sigEdges$geneC[[2]]))) == 0)
  expect_true(length(setdiff(rownames(subnetworksFull$geneC[[3]]), rownames(sigEdges$geneC[[3]]))) == 0)
})
test_that("[UNAGI] BuildSubnetworkU() function yields expected results",{
  
  # Construct a starting network, which will be modified.
  # Here, we expect genes A and B to be connected after 1 hop via TF2, genes A and D
  # to be connected after 1 hop via TF3, and genes A and C to be connected after 2
  # hops via TF4.
  startingNetwork <- data.frame(source = c(rep(c("gene1", "gene2", "gene3", "gene4"), 4)),
                                target = c(rep("geneA", 4), rep("geneB", 4), rep("geneC", 4), rep("geneD", 4)),
                                score = c(-3, 3, 5, -5, 0, 4, 0.0005, -0.5, 0, 0.0005, -1, 4, 0.0005, -2, 5, 3))
  addedNoise1 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  addedNoise2 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 10000)
  addedNoise3 <- data.frame(source = rep(1:100, 4), target = c(rep("geneA", 100), 
                                                         rep("geneB", 100),
                                                         rep("geneC", 100),
                                                         rep("geneD", 100)),
                            score = stats::rnorm(400) / 1000)
  fullNetwork1 <- rbind(startingNetwork, addedNoise1)
  rownames(fullNetwork1) <- paste(fullNetwork1$source, fullNetwork1$target, sep = "__")
  
  # Create a set of networks close to the original.
  fullNetwork2 <- fullNetwork1
  fullNetwork3 <- fullNetwork1
  fullNetwork2$score <- c(fullNetwork1$score - 0.00005)
  fullNetwork3$score <- c(fullNetwork1$score + 0.00005)
  
  # Paste together the networks.
  combinedNetwork <- fullNetwork1
  combinedNetwork[,4] <- fullNetwork2$score
  combinedNetwork[,5] <- fullNetwork3$score
  
  # Null distribution is a narrower version of the normal distribution.
  null <- rnorm(n = nrow(combinedNetwork) * 3) / 100
  
  # Verify that FindConnectionsForAllHopCountsU() and FindSignificantEdgesforhopU() are
  # working in tandem.
  # Obtain the overlaps for 1, 2, and 3 hops.
  expect_true(length(setdiff(c("gene2__geneA", "gene2__geneB", "gene3__geneA", "gene3__geneD", "gene4__geneC", "gene4__geneD"),
                             rownames(BuildSubnetworkU(geneSet = c("geneA", "geneB", "geneC", "geneD"),
                                        networks = list(fullNetwork1, fullNetwork2, fullNetwork3), 
                                        alpha = 0.1, hopConstraint = 4, nullDistribution = null)))) == 0)
  
  # Obtain the overlaps for only the first two hops, for only A and C.
  expect_true(length(setdiff(c("gene3__geneA", "gene3__geneD", "gene4__geneC", "gene4__geneD"),
                             rownames(BuildSubnetworkU(geneSet = c("geneA", "geneC"),
                                                      networks = list(fullNetwork1, fullNetwork2, fullNetwork3), 
                                                      alpha = 0.1, hopConstraint = 4, nullDistribution = null)))) == 0)
})
test_that("[UNAGI] RunUNAGI() function yields expected results",{
  
  # Check errors.
  expect_error(RunUNAGI(geneSet = 4, networks = list(), alpha = 0.5, hopConstraint = 5, nullDistribution = c(0,0,0)),
               paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
               "alpha and hopConstraint must be scalar numeric values."))
  expect_error(RunUNAGI(geneSet = "g1", networks = 67, alpha = 0.5, hopConstraint = 5, nullDistribution = c(0,0,0)),
               paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
                     "alpha and hopConstraint must be scalar numeric values."))
  expect_error(RunUNAGI(geneSet = "g1", networks = list(), alpha = "hi", hopConstraint = 5, nullDistribution = c(0,0,0)),
               paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
                     "alpha and hopConstraint must be scalar numeric values."))
  expect_error(RunUNAGI(geneSet = "g1", networks = list(), alpha = 0.5, hopConstraint = "hi", nullDistribution = c(0,0,0)),
               paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
                     "alpha and hopConstraint must be scalar numeric values."))
  expect_error(RunUNAGI(geneSet = "g1", networks = list(1,2,3), alpha = 0.5, hopConstraint = 5, nullDistribution = c(0,0,0)),
               "Each network must be a data frame.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA, blah = NA)), 
                           alpha = 0.5, hopConstraint = 5, nullDistribution = c(0,0,0)),
               paste("Each network must have transcription factors in the first column,",
               "target genes in the second column, and scores in the third column."))
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = -1, hopConstraint = 5, nullDistribution = c(0,0,0)), "alpha must be between 0 and 1, not including 0.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = 0, hopConstraint = 5, nullDistribution = c(0,0,0)), "alpha must be between 0 and 1, not including 0.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = 1.2, hopConstraint = 5, nullDistribution = c(0,0,0)), "alpha must be between 0 and 1, not including 0.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = 1, hopConstraint = -4, nullDistribution = c(0,0,0)), "hopConstraint must be an even number of at least 2.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = 1, hopConstraint = 7, nullDistribution = c(0,0,0)), "hopConstraint must be an even number of at least 2.")
  expect_error(RunUNAGI(geneSet = "g1", networks = list(data.frame(source = NA, target = NA, score = NA)), 
                           alpha = 1, hopConstraint = 4, nullDistribution = "hi"), "nullDistribution must be numeric.")
})