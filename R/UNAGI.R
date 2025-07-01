#' Given a set of genes of interest, full unipartite networks with scores (one network for each sample), a significance
#' cutoff for statistical testing, and a hop constraint, UNAGI finds a subnetwork of
#' significant edges connecting the genes.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' @param pValueChunks The number of chunks to split when calculating the p-value. This
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' @param pValueFile The file where the p-values should be saved. If NULL, they are not
#' saved and need to be recalculated.
#' @param loadPValues Whether p-values should be loaded from pValueFile or re-generated.
#' Default is FALSE.
#' @returns A unipartite subnetwork in the same format as the original networks.
#' @export

RunUNAGI <- function(nodeSet, network, alpha, hopConstraint, nullDistribution, 
                     verbose = FALSE, topX=NULL doFDRAdjustment = TRUE,
                     pValueChunks = 100, loadPValues = FALSE, pValueFile = "pvalues.RDS") {

  
  # === SAME AS BLOBFISH ===
 # Check for invalid inputs.
  if (!is.character(nodeSet) || !is.data.frame(network) || !is.numeric(alpha)) {
    stop("nodeSet must be character, network must be data frame, alpha must be numeric")
  } else if (!all(c("node1", "node2", "score") %in% colnames(network))) {
    stop("Network must have columns: node1, node2, score")
  } else if (alpha > 1 || alpha <= 0) {
    stop("alpha must be between 0 and 1, not including 0")
  } else if (!is.numeric(nullDistribution)) {
    stop("nullDistribution must be numeric")
  }
  
  # === SAME AS BLOBFISH ===
  # Calculate p-values (or load if requested)
  if (loadPValues) {
    pValues <- readRDS(pValueFile)
  } else {
    pValues <- CalculatePValuesUnipartite(network, nullDistribution, pValueChunks, doFDRAdjustment, pValueFile, verbose)
  }
  
  # Identify significant edges
  sigEdges <- network[pValues < alpha, ]
  sigEdges$pValue <- pValues[pValues < alpha]
  
  if (verbose) {
    message(paste("Retained", nrow(sigEdges), "significant edges"))
  }
  
  # Build clusters using greedy growth
  clusters <- BuildUnipartiteClusters(sigEdges, nodeSet, verbose)
  
  return(list(clusters = clusters, edges = sigEdges))
}
CalculatePValuesUnipartite <- function(network, nullDistribution, pValueChunks = 100,
                                       doFDRAdjustment = TRUE, pValueFile = "pvalues.RDS",
                                       verbose = FALSE) {

  
# === SAME AS BLOBFISH ===
  # Initialize p-values.
  pValues <- rep(NA, nrow(network))
  
  # Set the initial start and end indices.
  startIndex <- 1
  chunkSize <- ceiling(nrow(network) / pValueChunks)
  for (chunk in 1:pValueChunks) {
    endIndex <- min(startIndex + chunkSize - 1, nrow(network))
    
    scores <- network$score[startIndex:endIndex]
    pValues[startIndex:endIndex] <- sapply(scores, function(score) {
      sum(nullDistribution >= score) / length(nullDistribution)
    })
    
    if (verbose) {
      message(paste("Processed chunk", chunk, "of", pValueChunks))
    }
    
    if (endIndex == nrow(network)) break
    startIndex <- endIndex + 1
  }
  
  if (doFDRAdjustment) {
    pValues <- p.adjust(pValues, method = "fdr")
  }
  
  names(pValues) <- paste(network$node1, network$node2, sep = "__")
  
  if (!is.null(pValueFile)) {
    saveRDS(pValues, pValueFile)
  }
  
  return(pValues)
}

BuildUnipartiteClusters <- function(sigEdges, nodeSet, verbose = FALSE) {
  # === Analogy to BLOBFISH BuildSubnetwork + FindConnections ===
  
  clusters <- list()
  clusterMap <- list()
  clusterID <- 1
  
  for (i in 1:nrow(sigEdges)) {
    n1 <- sigEdges$node1[i]
    n2 <- sigEdges$node2[i]
    
    in1 <- n1 %in% names(clusterMap)
    in2 <- n2 %in% names(clusterMap)
    
    if (!in1 && !in2) {
      # New cluster
      clusters[[clusterID]] <- c(n1, n2)
      clusterMap[[n1]] <- clusterID
      clusterMap[[n2]] <- clusterID
      clusterID <- clusterID + 1
    } else if (in1 && !in2) {
      cid <- clusterMap[[n1]]
      clusters[[cid]] <- unique(c(clusters[[cid]], n2))
      clusterMap[[n2]] <- cid
    } else if (!in1 && in2) {
      cid <- clusterMap[[n2]]
      clusters[[cid]] <- unique(c(clusters[[cid]], n1))
      clusterMap[[n1]] <- cid
    } else {
      # Merge clusters
      cid1 <- clusterMap[[n1]]
      cid2 <- clusterMap[[n2]]
      if (cid1 != cid2) {
        clusters[[cid1]] <- unique(c(clusters[[cid1]], clusters[[cid2]]))
        for (node in clusters[[cid2]]) {
          clusterMap[[node]] <- cid1
        }
        clusters[[cid2]] <- NULL
      }
    }
  }
  
  # Filter to include clusters with at least one node of interest
  clusters <- Filter(function(cl) any(cl %in% nodeSet), clusters)
  
  if (verbose) {
    message(paste("Final number of clusters with nodes of interest:", length(clusters)))
  }
  
  return(clusters)
}


