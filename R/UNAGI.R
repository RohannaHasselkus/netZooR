#' Given a set of genes of interest, full unipartite networks with scores (one network for each sample), a significance
#' cutoff for statistical testing, and a hop constraint, UNAGI finds a subnetwork of
#' significant edges connecting the genes.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' saved and need to be recalculated.
#' Default is FALSE.
#' @returns A unipartite subnetwork in the same format as the original networks.
#' @export
RunUNAGI <- function(nodeSet, network, alpha, hopConstraint, 
                     verbose = FALSE, topX=NULL) {
  
  # === SAME AS BLOBFISH ===
 # Check for invalid inputs.
  if (!is.character(nodeSet) || !is.data.frame(network) || !is.numeric(alpha)) {
    stop("nodeSet must be character, network must be data frame, alpha must be numeric")
  } else if (!all(c("node1", "node2", "score") %in% colnames(network))) {
    stop("Network must have columns: node1, node2, score")
  } else if (alpha > 1 || alpha <= 0) {
    stop("alpha must be between 0 and 1, not including 0")
  }
    
  # Build clusters using greedy growth
  clusters <- BuildUnipartiteClusters(sigEdges, nodeSet, verbose)
  
  return(list(clusters = clusters, edges = sigEdges))
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
      clusterMap[[n2]] <- cid    } else if (!in1 && in2) {
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
#THE REST IS THE EXACT SAME AS BLOBFISH
#' Find the subnetwork of significant edges n / 2 hops away from each gene.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param combinedNetwork A concatenation of n PANDA-like networks with the following format:
#' @param pValues The p-values for all edges.
#' @param hopConstraint The maximum number of hops to be considered for a gene.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
FindSignificantEdgesForHop <- function(geneSet, combinedNetwork, hopConstraint, pValues,
                                       verbose = FALSE, topX = NULL){
  # Build the significant subnetwork for each gene, up to the hop constraint.
  uniqueGeneSet <- sort(unique(geneSet))
  geneSubnetworks <- lapply(uniqueGeneSet, function(gene){
    
    # Get all significant edges for a 1-hop subnetwork.
    if(verbose == TRUE){
      message(paste("Evaluating hop 1 for gene", gene))
    }
    subnetwork1Hop <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                    pValues = pValues, 
                                                    startingNodes = gene,
                                                    nodesToExclude = c(),
                                                    verbose = verbose,
                                                    topX = topX)
    
    # Set the starting and excluded set for the next hop.
    startingNodes <- unique(subnetwork1Hop$)
    topXNew <- NULL
    if(!is.null(topX)){
      topXNew <- topX * length(startingNodes <- unique(c(subnetwork1Hop[,1], subnetwork1Hop[,2])))
    }
    excludedSubset <- gene
    
    # Add to the list of all subnetworks.
    allSubnetworksForGene <- list(subnetwork1Hop)
    
    # Loop until we reach the maximum number of hops or there are no new edges
    # to traverse.
    hop <- 2
    while(hop <= hopConstraint && length(startingNodes) > 0){
      
      # If we are on an even hop, start from transcription factors.
      # If we are on an odd hop, start from genes.
      if(hop %% 2 == 0){
        
        # Find all significant edges in the next hop.
        if(verbose == TRUE){
          message(paste("Evaluating hop", hop, "for gene", gene))
        }
        subnetworkHops <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                        pValues = pValues, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset,
                                                     
                                                        verbose = verbose,
                                                        topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(subnetworkHops$gene), excludedSubset)
        if(!is.null(topX)){
          topXNew <- topX * length(startingNodes)
        }
      }else{
        
        # Find all significant edges in the next hop.
        if(verbose == TRUE){
          message(paste("Evaluating hop", hop, "for gene", gene))
        }
        subnetworkHops <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                        pValues = pValues, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset, 
                                                        verbose = verbose, 
                                                        topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(subnetworkHops<- unique(c(subnetwork1Hop[,1], subnetwork1Hop[,2]))), excludedSubset)
        if(!is.null(topX)){
          topXNew <- topX * length(startingNodes)
        }
      }
      
      # Add to the list.
      allSubnetworksForGene[[length(allSubnetworksForGene) + 1]] <- subnetworkHops
      
      # Increment hops.
      hop <- hop + 1
    }
    return(allSubnetworksForGene)
  })
  
  # Add the names of the genes.
  names(geneSubnetworks) <- uniqueGeneSet
  return(geneSubnetworks)
}
#' Find all significant edges adjacent to the starting nodes, excluding the nodes
#' specified.
#' @param networks A concatenation of n PANDA-like networks with the following format:
#' @param pValues The p-values from the original network.
#' @param startingNodes The list of nodes from which to start.
#' @param nodesToExclude The list of nodes to exclude from the search.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
SignificantBreadthFirstSearch <- function(networks, pValues, startingNodes,
                                          nodesToExclude, 
                                          verbose = FALSE, topX = NULL){
  # Check that provided nodes overlap with the networks.
  if((length(setdiff(startingNodes, c(networks[,1],networks[,2])) > 0)){
  
  # Identify genes and transcription factors to test, based on which of these we are
  # starting from.
 genesToTest <- setdiff(unique(c(networks[,1],networks[,2]), nodesToExclude)
  
  # Construct all edges to test based on the combination of these.
  srcGeneLongList <- rep(genesToTest, length(genesToTest))
tgtGeneLongList <- unlist(lapply(genesToTest, function(gene){
  return(rep(gene, length(genesToTest)))
}))
subnetwork <- networks[
  (networks$node1 %in% startingNodes | networks$node2 %in% startingNodes) &
  !(networks$node1 %in% nodesToExclude | networks$node2 %in% nodesToExclude),
, ]
  
  # For each edge, measure its significance.
  subnetwork <- networks
  if(length(allEdges) > 0){
    
    # If topX is specified, filter again.
    significantEdges <- allEdges
    if(!is.null(topX) && length(allEdges) > topX){
      whichTopX <- order(pValues[allEdges])[1:topX]
      significantEdges <- allEdges[whichTopX]
    }
    
    # Return the edges meeting alpha.
    subnetwork <- networks[significantEdges, c(1:2)]
    if(verbose == TRUE){
      message(paste("Retained", length(significantEdges), "edges"))
    }
  }
  
  # Return the subnetwork.
  return(subnetwork)
}
#' For all hop counts up to the maximum, find subnetworks connecting each pair of
#' genes by exactly that number of hops. For instance, find each 
#' containing only the significant edges meeting the hop count criteria and
#' where each network is a data frame with the following format:
#' @param verbose Whether or not to print detailed information about the run.
FindConnectionsForAllHopCounts <- function(subnetworks, verbose = FALSE){
  
  # Find a subnetwork for each hop count.
  hopCountSubnetworks <- lapply(1:length(subnetworks[[1]]), function(hops){
    
    # For each pair of genes, find the subnetworks for this number of hops.
    geneSpecificHopCountSubnetwork <- lapply(1:(length(names(subnetworks))-1), function(i){
      genePairSpecificHopCountSubnetwork <- lapply((i+1):length(names(subnetworks)), function(j){
        
        # Get the subnetworks for the number of hops of interest.
        gene1 <- names(subnetworks)[i]
        gene2 <- names(subnetworks)[j]
        connectingSubnetwork <- data.frame(source = NA, target = NA)[0,]
        
        # If there were no edges at this hop count for one or both genes,
        # do not evaluate.
        if(length(subnetworks[[gene1]]) >= hops && length(subnetworks[[gene2]]) >= hops){
          subnetwork1 <- subnetworks[[gene1]][[hops]]
          subnetwork2 <- subnetworks[[gene2]][[hops]]
          
          # Initialize overlapping subnetwork.
          sourceToRecurse1 <- c()
          sourceToRecurse2 <- c()
          targetToRecurse1 <- c()
          targetToRecurse2 <- c()
          
          # If the number of hops is even, add edges from genes that overlap
          # If the number of hops is odd, add edges from transcription factors that overlap.
          if(hops %% 2 == 0){
            overlappingGenes <- intersect(subnetwork1$gene, subnetwork2$gene)
            if(verbose == TRUE){
              message(paste("Hop", hops, "-", length(overlappingGenes), "overlapped between", gene1, "and", gene2))
            }
            whichSubnet1Gene <- which(subnetwork1$gene %in% overlappingGenes)
            whichSubnet2Gene <- which(subnetwork2$gene %in% overlappingGenes)
            tfsToRecurse1 <- unique(subnetwork1[whichSubnet1Gene, "tf"])
            tfsToRecurse2 <- unique(subnetwork2[whichSubnet2Gene, "tf"])
            connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichSubnet1Gene,],
                                          subnetwork2[whichSubnet2Gene,])
          }else{
            overlappingTF <- intersect(subnetwork1$tf, subnetwork2$tf)
            if(verbose == TRUE){
              message(paste("Hop", hops, "-", length(overlappingTF), "overlapped between", gene1, "and", gene2))
            }
            whichSubnet1TF <- which(subnetwork1$tf %in% overlappingTF)
            whichSubnet2TF <- which(subnetwork2$tf %in% overlappingTF)
            genesToRecurse1 <- unique(subnetwork1[whichSubnet1TF, [,2]])
            genesToRecurse2 <- unique(subnetwork2[whichSubnet2TF, [,2]])
            connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichSubnet1TF,],
                                          subnetwork2[whichSubnet2TF,])
          }
          
          # Recurse back over the number of hops.
          if(hops-1 >= 1){
            for(hop in (hops-1):1){
              subnetwork1 <- subnetworks[[gene1]][[hop]]
              subnetwork2 <- subnetworks[[gene2]][[hop]]
              
              # If the current number of hops is even, add edges from TFs connected to genes of interest.
              # If the current number of hops is odd, add edges from genes connected to TFs of interest.
              if(hop %% 2 == 0){
                whichTFConnectedToGene1 <- which(subnetwork1$gene %in% genesToRecurse1)
                whichTFConnectedToGene2 <- which(subnetwork2$gene %in% genesToRecurse2)
                tfsToRecurse1 <- unique(subnetwork1[whichTFConnectedToGene1, "tf"])
                tfsToRecurse2 <- unique(subnetwork2[whichTFConnectedToGene2, "tf"])
                connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichTFConnectedToGene1,],
                                              subnetwork2[whichTFConnectedToGene2,])
              }            }
          }
        }
        
        # Return the subnetwork, which should now contain all of the edges connecting the
        # gene pair at the prespecified number of hops.
        return(connectingSubnetwork)
      })
      # Bind together the subnetwork for each gene pair.
      connectingSubnetworkAll <- do.call(rbind, genePairSpecificHopCountSubnetwork)
      return(connectingSubnetworkAll)
    })
    
    # Bind together the subnetworks for each gene.
    return(do.call(rbind, geneSpecificHopCountSubnetwork))
  })
  
  # Bind together the subnetworks for each hop count.
  compositeSubnetwork <- do.call(rbind, hopCountSubnetworks)
  compositeSubnetworkEdges <- paste(compositeSubnetwork$node1, compositeSubnetwork$node2, sep = "__")
  uniqueEdges <- sort(unique(compositeSubnetworkEdges))
  compositeSubnetworkDedup <- do.call(rbind, lapply(uniqueEdges, function(edge){
    whichFirstEdge <- which(compositeSubnetworkEdges == edge)[1]
    return(compositeSubnetwork[whichFirstEdge,])
  }))
  rownames(compositeSubnetworkDedup) <- uniqueEdges
  
  # Remove all genes connected to a single transcription factor. These genes were
  # added because they are regulated by a transcription factor that co-regulates
  # two seed genes. Similarly, remove all transcription factors connected to a 
  # single gene.
  geneCounts <- table(compositeSubnetworkDedup$gene)
  tfCounts <- table(compositeSubnetworkDedup$tf)
  genesToRemove <- names(geneCounts)[which(geneCounts == 1)]
  genesToRemove <- setdiff(genesToRemove, names(subnetworks))
  tfsToRemove <- names(tfCounts)[which(tfCounts == 1)]
  compositeSubnetworkDedup <- compositeSubnetworkDedup[which(compositeSubnetworkDedup$gene %in% setdiff(names(geneCounts), 
                                                                                                        genesToRemove)),]
  compositeSubnetworkDedup <- compositeSubnetworkDedup[which(compositeSubnetworkDedup$tf %in% setdiff(names(tfCounts), 
                                                                                                      tfsToRemove)),]
  return(compositeSubnetworkDedup)
}
#' Plot the networks, using different colors for transcription factors, genes of interest,
#' and additional genes.
#' @param network A data frame with the following format:
#' tf,gene
#' @param genesOfInterest Which genes of interest to highlight
#' @param tfColor Color for the transcription factors
#' @param geneColorMapping Color mapping from a set of genes to a color. The
#' nodes and edges connected to them will be this color. If NULL, all genes and
#' their edges will be gray. The format is a data frame, where the first column ("gene")
#' is the name of the gene and the second ("color") is the color.
#' @param nodeSize Size of node
#' @param edgeWidth Width of edges
#' @param vertexLabels Which vertex labels to include. By default, none are included.
#' @param vertexLabelSize The size of label to use for the vertex, as a fraction of the default.
#' @param vertexLabelOffset Number of pixels in the offset when plotting labels.
#' Default is TRUE.
#' @export
PlotNetwork <- function(network, genesOfInterest,
                        tfColor = "blue", nodeSize = 1,
                        edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.7,
                        vertexLabelOffset = 0.5, geneColorMapping = NULL){
  
  # Convert from factor to character.
  network$node1 <- as.character(network$node1
  network$node2 <- as.character(network$node2)
  
  # Set the node attributes.
  uniqueNodes <- unique(c(network$node1, network$node2))
  nodeAttrs <- data.frame(node = uniqueNodes,
                          color = rep("gray", length(uniqueNodes)),
                          size = rep(nodeSize, length(uniqueNodes)),
                          frame.width = rep(0, length(uniqueNodes)),
                          label.color = "black", label.cex = vertexLabelSize,
                          label.dist = vertexLabelOffset)
  rownames(nodeAttrs) <- uniqueNodes
  
  # Add colors.
  nodeAttrs$color <- "gray"
if (!is.null(geneColorMapping)) {
  rownames(geneColorMapping) <- geneColorMapping$gene
  nodeAttrs[rownames(nodeAttrs) %in% rownames(geneColorMapping), "color"] <- 
    geneColorMapping[rownames(nodeAttrs)[rownames(nodeAttrs) %in% rownames(geneColorMapping)], "color"]
}

  # Add gene colors.
  rownames(geneColorMapping) <- geneColorMapping$gene
  geneColorMapping <- geneColorMapping[intersect(rownames(geneColorMapping), uniqueNodes),]
  if(!is.null(geneColorMapping)){
    nodeAttrs[rownames(geneColorMapping), "color"] <- geneColorMapping$color
  }
  # Add edge attributes.
  if(!is.null(geneColorMapping)){
    for(gene in rownames(geneColorMapping)){
      network[which(network$node2 == [,2]), "color"] <- geneColorMapping[gene, "color"]
    }
  }
  network$width <- edgeWidth
  # Create a graph object.
  graph <- igraph::graph_from_data_frame(network, vertices = nodeAttrs, directed = FALSE)
  V(graph)$type <- V(graph)$name %in% network$node1
  
  # Plot.
  labels <- V(graph)$name
  whichEmpty <- which(labels %in% setdiff(labels, vertexLabels))
  labels[whichEmpty] <- rep(NA, length(whichEmpty))
 }

