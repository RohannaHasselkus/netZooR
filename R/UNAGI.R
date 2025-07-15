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
  #this entire body was edited Jul 13
  if (!is.character(nodeSet) || !is.data.frame(network) || !is.numeric(alpha))
    stop("Wrong input type! geneSet must be a character vector. networks must be a list.",
         "alpha and hopConstraint must be scalar numeric values.")
  if (!all(c("source", "target", "score") %in% colnames(network)))
    stop("Each network must have transcription factors in the first column,",
         "target genes in the second column, and scores in the third column.")
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be between 0 and 1, not including 0")
  
  #significance filter
  sigEdges <- network[, c("source", "target")]
  rownames(sigEdges) <- paste(sigEdges$source, sigEdges$target, sep = "__")
  
  ##Build clusters exactly like BLOBFISH 
  clusters <- BuildUnipartiteClusters(sigEdges, nodeSet, verbose)
  
  list(clusters = clusters, edges = sigEdges)
}

BuildUnipartiteClusters <- function(sigEdges, nodeSet, verbose = FALSE) {
  # === Analogy to BLOBFISH BuildSubnetwork + FindConnections ===
  
  clusters <- list()
  clusterMap <- list()
  clusterID <- 1
  
  for (i in 1:nrow(sigEdges)) {
    n1 <- sigEdges$source[i]
    n2 <- sigEdges$target[i]
    
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
      clusterMap[[n2]] <- cid} else if (!in1 && in2) {
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

#' Find the subnetwork of significant edges connecting the genes.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param networks A list of  unipartite (PANDA-like) networks, where each network is a data frame with the following format:
#' tf,gene,score
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' Must be an even number.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' saved and need to be recalculated.
#' Default is FALSE.
#' @returns A  unipartite subnetwork in the same format as the original networks.
BuildSubnetworkU <- function(geneSet, networks, alpha, hopConstraint, nullDistribution,
                            verbose = FALSE, topX = NULL, doFDRAdjustment = TRUE){
  
  # Name edges for each network.
  combinedNetwork <- networks
  if(!is.data.frame(combinedNetwork)){
    networksNamed <- lapply(networks, function(network){
      rownames(network) <- paste(network[,2], network[,1], sep = "__")
      return(network)
    })
    # Paste together the networks.
    combinedNetwork <- networksNamed[[1]]
    for(i in 2:length(networksNamed)){
      combinedNetwork[,2+i] <- networksNamed[[i]]$score
    }
  }
 
  
  # Compute significance mask (placeholder)
  whichSig <- rep(TRUE, nrow(combinedNetwork))
   
  # Subset the network.
  significantEdges <- rownames(combinedNetwork)[whichSig]
  subnetwork <- combinedNetwork[significantEdges, c(1:2)]
  genesWithNoSigEdges <- setdiff(geneSet, subnetwork[,1])
  geneSet <- intersect(geneSet, subnetwork[,1])
  if(verbose == TRUE){
    message(paste("The following genes had no significant edges:", paste(genesWithNoSigEdges, collapse = ",")))
    message(paste("Retained", length(significantEdges), "out of", length(rownames(combinedNetwork)), "edges"))
  }
  
  # For each gene, find the significant edges from each hop.
  significantSubnetworks <- FindSignificantEdgesForHop(geneSet = geneSet,
                                                       combinedNetwork = subnetwork,
                                                       hopConstraint = hopConstraint / 2,
                                                       verbose = verbose, topX = topX)
}

#anlogous to BLOBFISH
#' Find the subnetwork of significant edges n / 2 hops away from each gene.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param combinedNetwork A concatenation of n PANDA-like networks with the following format:
#' @param hopConstraint The maximum number of hops to be considered for a gene.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
FindSignificantEdgesforhopU <- function(geneSet, combinedNetwork, hopConstraint,
                                        verbose = FALSE, topX = NULL){
  # Build the significant subnetwork for each gene, up to the hop constraint.
  uniqueGeneSet <- sort(unique(geneSet))
  geneSubnetworks <- lapply(uniqueGeneSet, function(gene){
    
    # Get all significant edges for a 1-hop subnetwork.
    if(verbose == TRUE){
      message(paste("Evaluating hop 1 for gene", tf))
    }
    subnetwork1Hop <- SignificantBreadthFirstSearchU(networks = combinedNetwork,
                                                     startingNodes = gene,
                                                     nodesToExclude = c(),
                                                     verbose = verbose,
                                                     topX = topX)
    
    # Set the starting and excluded set for the next hop.
    startingNodes <- unique(subnetwork1Hop)
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
        subnetworkHops <- SignificantBreadthFirstSearchU(networks = combinedNetwork,
                                                         startingNodes = startingNodes,
                                                         nodesToExclude = excludedSubset,
                                                         verbose = verbose,
                                                         topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(subnetworkHops[,1]), excludedSubset)
        if(!is.null(topX)){
          topXNew <- topX * length(startingNodes)
        }
      }else{
        
        # Find all significant edges in the next hop.
        if(verbose == TRUE){
          message(paste("Evaluating hop", hop, "for gene", gene))
        }
        subnetworkHops <- SignificantBreadthFirstSearchU(networks = combinedNetwork,
                                                         startingNodes = startingNodes,
                                                         nodesToExclude = excludedSubset,
                                                         verbose = verbose,
                                                         topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(c(subnetworkHops[,1], subnetworkHops[,2])), excludedSubset)
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
#' @param networks A PANDA-like network
#' @param startingNodes The list of nodes from which to start.
#' @param nodesToExclude The list of nodes to exclude from the search.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
SignificantBreadthFirstSearchU <- function(networks, startingNodes,
                                           nodesToExclude,
                                           verbose = FALSE, topX = NULL){

  # Check that provided nodes overlap with the networks.
  if(length(setdiff(startingNodes, c(networks[,1], networks[,2]))) > 0){
    stop("ERROR: Starting nodes do not overlap with network nodes")
  }
  if(length(setdiff(nodesToExclude, c(networks[,1], networks[,2]))) > 0){
    stop("ERROR: List of nodes to exclude does not overlap with network nodes")
  }
  if(length(intersect(startingNodes, nodesToExclude)) > 0){
    stop("ERROR: Starting nodes cannot overlap with nodes to exclude")
  }
  
  
  # Identify genes and transcription factors to test, based on which of these we are
  # starting from.
  genesToTest <- setdiff(unique(c(networks[,1],networks[,2])), nodesToExclude)
  
  # Construct all edges to test based on the combination of these.
  srcGeneLongList <- rep(genesToTest, length(genesToTest))
  tgtGeneLongList <- unlist(lapply(genesToTest, function(gene){
    return(rep(gene, length(genesToTest)))
  }))
  subnetwork <- networks[
    (networks$source %in% startingNodes | networks$target %in% startingNodes) &
      !(networks$source %in% nodesToExclude | networks$target %in% nodesToExclude),
    , ]
  
  # For each edge, measure its significance.
  subnetwork <- networks
  allEdges <- rownames(networks)
  if(length(allEdges) > 0){
    
    # If topX is specified, filter again.
    significantEdges <- allEdges
    if(!is.null(topX) && length(allEdges) > topX){
      whichTopX <- order(allEdges)[1:topX]
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
FindConnectionsForAllHopCountsU <- function(subnetworks, verbose = FALSE){
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
          
          #alreadyExploredGenes <- c(gene1, gene2)

          # Initialize overlapping subnetwork.
          geneToRecurse1 <- c()
          geneToRecurse2 <- c()
          
          # Add edges from genes that overlap.
          overlappingGenes <- intersect(subnetwork1[, 2], subnetwork2[, 2])
          if(verbose == TRUE){
            message(paste("Hop", hops, "-", length(overlappingGenes), "overlapped between", gene1, "and", gene2))
          }
          whichSubnet1Gene <- which(subnetwork1[,2] %in% overlappingGenes)
          whichSubnet2Gene <- which(subnetwork2[,2] %in% overlappingGenes)
          geneToRecurse1 <- unique(subnetwork1[whichSubnet1Gene, 1])
          geneToRecurse2 <- unique(subnetwork2[whichSubnet2Gene, 1])
          if(hops == 3 && "gene3" %in% overlappingGenes){
            print(overlappingGenes)
            print(subnetwork1[whichSubnet1Gene,])
            print(subnetwork2[whichSubnet2Gene,])
            print(geneToRecurse1)
            print(geneToRecurse2)
          }
          
          # Subset genes to recurse
          #geneToRecurse1 <- setdiff(geneToRecurse1, alreadyExploredGenes)
          #geneToRecurse2 <- setdiff(geneToRecurse2, alreadyExploredGenes)
          # Update explored genes
          #alreadyExploredGenes <- union(alreadyExploredGenes, c(geneToRecurse1, geneToRecurse2))
          
          # Add overlapping edges to subnetwork
          connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichSubnet1Gene,],
                                        subnetwork2[whichSubnet2Gene,])
          
          # Recurse back over the number of hops.
          if(hops-1 >= 1){
            for(hop in (hops-1):1){
              subnetwork1 <- subnetworks[[gene1]][[hop]]
              subnetwork2 <- subnetworks[[gene2]][[hop]]
              
              # If the current number of hops is even, add edges from genes connected to genes of interest.
              whichGeneConnectedToGene1 <- which(subnetwork1[,2] %in% geneToRecurse1)
              whichGeneConnectedToGene2 <- which(subnetwork2[,2] %in% geneToRecurse2)
              
              genesToRecurse1 <- unique(subnetwork1[whichGeneConnectedToGene1, 1])
              genesToRecurse2 <-unique(subnetwork2[whichGeneConnectedToGene2, 1])
              if(hops == 3 && "gene3" %in% overlappingGenes){
                print(hop)
                print(overlappingGenes)
                print(subnetwork1[whichGeneConnectedToGene1,])
                print(subnetwork2[whichGeneConnectedToGene2,])
                print(subnetwork1)
                print(subnetwork2)
                print(geneToRecurse1)
                print(geneToRecurse2)
              }
              #genesToRecurse1 <- setdiff(unique(subnetwork1[whichGeneConnectedToGene1, 1]), alreadyExploredGenes)
              #genesToRecurse2 <- setdiff(unique(subnetwork2[whichGeneConnectedToGene2, 1]), alreadyExploredGenes)
              #alreadyExploredGenes <- union(alreadyExploredGenes, c(genesToRecurse1, genesToRecurse2))
              
              connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichGeneConnectedToGene1,],
                                            subnetwork2[whichGeneConnectedToGene2,])
            }
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
  
  #Fix ordering nd filtering
  compositeSubnetwork <- do.call(rbind, hopCountSubnetworks)
  colnames(compositeSubnetwork) <- c("source", "target")

  # remove duplicates in real time
  edgeKeys <- paste(compositeSubnetwork$source,
                    compositeSubnetwork$target, sep = "__")
  edgeKeysReverse <- paste(compositeSubnetwork$target,
                           compositeSubnetwork$source, sep = "__")
  compositeSubnetworkDedup <- compositeSubnetwork
  
  # Remove all duplicated edges.
  i = 1
  len <- length(edgeKeys)
  while(i < len){
    whichToRemove <- setdiff(which(edgeKeys == edgeKeys[i]), i)
    if(length(whichToRemove) > 0){
      compositeSubnetworkDedup <- compositeSubnetworkDedup[-whichToRemove,]
      edgeKeys <- edgeKeys[-whichToRemove]
      edgeKeysReverse <- edgeKeysReverse[-whichToRemove]
      len <- length(edgeKeys)
    }
    i <- i + 1
  }

  # Remove all reversed edges.
  i = 1
  len <- length(edgeKeys)
  while(i < len){
    if(edgeKeysReverse[i] %in% edgeKeys){
      whichToRemove <- which(edgeKeys == edgeKeysReverse[i])
      compositeSubnetworkDedup <- compositeSubnetworkDedup[-whichToRemove,]
      edgeKeys <- edgeKeys[-whichToRemove]
      edgeKeysReverse <- edgeKeysReverse[-whichToRemove]
      len <- length(edgeKeys)
    }
    i <- i + 1
  }
  
  # Set row names.
  rownames(compositeSubnetworkDedup) <- paste(compositeSubnetworkDedup[,1], compositeSubnetworkDedup[,2], sep = "__")
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
PlotNetworkU <- function(network, genesOfInterest,
                         tfColor = "blue", nodeSize = 1,
                         edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.7,
                         vertexLabelOffset = 0.5, geneColorMapping = NULL){
  
  # Convert from factor to character.
  network$source <- as.character(network$source)
  network$target <- as.character(network$target)
  
  # Set the node attributes.
  uniqueNodes <- unique(c(network$source, network$target))
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
    rownames(geneColorMapping) <- geneColorMapping[,1]
    nodeAttrs[rownames(nodeAttrs) %in% rownames(geneColorMapping), "color"] <-
      geneColorMapping[rownames(nodeAttrs)[rownames(nodeAttrs) %in% rownames(geneColorMapping)], "color"]
  }
  
  # Add gene colors.
  rownames(geneColorMapping) <- geneColorMapping[,1]
  geneColorMapping <- geneColorMapping[intersect(rownames(geneColorMapping), uniqueNodes),]
  if(!is.null(geneColorMapping)){
    nodeAttrs[rownames(geneColorMapping), "color"] <- geneColorMapping$color
  }
  # Add edge attributes.
  if(!is.null(geneColorMapping)){
    for(gene in rownames(geneColorMapping)){
      network[which(network$target == 2), "color"] <- geneColorMapping[gene, "color"]
    }
  }
  network$width <- edgeWidth
  # Create a graph object.
  graph <- igraph::graph_from_data_frame(network, vertices = nodeAttrs, directed = FALSE)
  V(graph)$type <- V(graph)$name %in% network$source
  
  # Plot.
  labels <- V(graph)$name
  whichEmpty <- which(labels %in% setdiff(labels, vertexLabels))
  labels[whichEmpty] <- rep(NA, length(whichEmpty))
}
