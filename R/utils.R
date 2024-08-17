# Nomalized data
centerCounts <- function(obj,
                         doInChunks=TRUE,
                         chunkSize=1000){
  if(!class(obj) %in% c("SummarizedExperiment","RangedSummarizedExperiment","dgCMatrix","dgeMatrix","Matrix"))
    stop("Supplied object must be either of class SummarizedExperiment or sparse Matrix ..\n")
  
  if(ncol(obj) > 10000)
    doInChunks <- TRUE
  
  if(doInChunks){
    cat("Centering counts for cells sequentially in groups of size ",
        chunkSize, " ..\n\n")
    starts <- seq(1,ncol(obj),chunkSize)
  } else{
    starts <- 1
  }
  
  counts.l <- list()
  
  for(i in 1:length(starts)){
    beginning <- starts[i]
    if(i==length(starts)) # If it's the last one
    {
      ending <- ncol(obj)
    } else {
      ending <- starts[i]+chunkSize-1
    }
    
    cat("Computing centered counts for cells: ",beginning," to ", ending,"..\n")
    
    if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
      m <- SummarizedExperiment::assay(obj[, beginning:ending])} else {
        m <- obj[,beginning:ending] # Assumes Matrix format
      }
    cellMeans <- Matrix::colMeans(m)
    cat("Computing centered counts per cell using mean reads in features ..\n\n")
    # Center cell counts based on its mean RIP count
    cCounts <- Matrix::t(Matrix::t(m)/cellMeans)
    
    counts.l[[i]] <- cCounts
    
    gc()
  }
  
  cat("Merging results..\n")
  centered.counts <- do.call("cbind",counts.l)
  cat("Done!\n")
  
  if(class(obj) == "RangedSummarizedExperiment" | class(obj)=="SummarizedExperiment"){
    SummarizedExperiment::assay(obj) <- centered.counts
    return(obj)
  } else {
    return(centered.counts)
  }
}

#Gene peak overlaps
genePeakOv <- function(ATAC.se, # SummarizedExperiment object of scATAC data
                       RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                       genome, # Must be one of "hg19", "mm10", "hg38", "macFas5"
                       geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                       windowPadSize=100000, # base pairs padded on either side of gene TSS
                       proPadSize = 2000 # base pairs padded on either side of gene TSS for enhancer
                        
) {
  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))
  
  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")
  
  message("Assuming paired scATAC/scRNA-seq data ..")
  
  ATACmat <- assay(ATAC.se) # Rownames preserved
  
  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")
  
  # Check for peaks/genes with 0 accessibility/expression
  
  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }
  
  
  peakRanges <- granges(ATAC.se) # Peak ranges
  names(peakRanges) <- rownames(ATAC.se)
  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }
  
  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")
  
  
  if (!genome %in% c("hg19", "hg38", "mm10", "macFas5"))
    stop("You must specify one of hg19, hg38, mm10, or macFas5 as a genome build for currently supported TSS annotations..\n")

  if(genome %in% c("hg19", "hg38", "mm10", "macFas5")){
    load(paste0('../data/', genome, '_refSeq.Rdata'))
    TSSg <- get(paste0(genome, 'TSSRanges'))
  }
  
  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)
  
  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")
    
    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }
    
    TSSg <- TSSg[geneList]
  }
  
  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")
  
  # Pad promoter by this much *either side*
  proflank <- GenomicRanges::flank(TSSg,
                                   width = proPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping promoter-gene pairs ..\n")
  geneProOv <- findOverlaps(query = proflank,subject = peakSummits)
  numpgPairs <- length(geneProOv)
  
  cat("Found ",numpgPairs,"total promoter-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene promoter window: ",length(unique(subjectHits(geneProOv))),"\n")
  cat("Number of gene promoter windows that overlap any peak summit: ",length(unique(queryHits(geneProOv))),"\n\n")
  peakSummits <- peakSummits[-unique(subjectHits(geneProOv))]
  
  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")
  
  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]
  
  
  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)
  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlapPairs(query = TSSflank,subject = peakSummits)
  genePeakOv <- data.frame(genePeakOv@first,genePeakOv@second)
  numPairs <- nrow(genePeakOv)
  
  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")
  
  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(genePeakOv$Peak)),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(genePeakOv$gene_name)),"\n\n")
  
  return(genePeakOv)

}

# Gene peak correlation
PeakGeneCor <- function(ATAC, # Normalized reads in peaks counts (rownames should  be "Peak1","Peak2" etc.)
                        RNA, # Normalized gene expression counts
                        OV, # Gene TSS - Peak overlap pairs object (Genes: query, Peaks: subject)
                        peakRanges,
                        seed = 2022,
                        mtd="spearman") {
  
  
  geneIndices <- as.character(OV$gene_name)
  peakIndices <- as.character(OV$Peak)
  
  uniquegenes <- unique(geneIndices)
  uniquepeaks <- unique(peakIndices)
  
  
  idy = uniquepeaks
  ovlpn = length(idy)
  idyexcl = which(!1:length(peakRanges) %in% idy)
  gene.name = uniquegenes
  
  gene.val = RNA[uniquegenes,];
  
  data.raw = ATAC
  peak.raw = peakRanges
  # use set
  data.use = ATAC[idy,,drop = FALSE];
  peak.use = peak.raw[uniquepeaks];
  
  corrFunc <- function(var1, var2, method) {
    result = cor.test(var1, var2, method = method)
    data.frame(result[c("estimate","p.value","statistic")], stringsAsFactors=FALSE)
  }
  
  #-----------------------------------
  # cat("Real set: performing statitical test on... \n", file = stderr())
  peaks.id = seq(nrow(data.use));
  corr = lapply(peaks.id, function(t){
    corrFunc(as.numeric(gene.val), as.numeric(data.use[t,]), mtd)
  }) %>% bind_rows()
  
  peak.use$estimate <- corr[, "estimate"]
  peak.use$statistic <- corr[, "statistic"]
  peak.use$method <- mtd
  peak.use$Pval <- corr[, "p.value"]
  peak.use$FDR <- p.adjust(peak.use$Pval, method = "BH")
  peak.use$class <- "corr"
  peak.use$Gene <- gene.name
  
  # cat("Done ... \n", file = stderr())
  
  #############################
  # OUTPUT
  return( as.data.frame(peak.use) );
}


# build enhancer network for each gene
Network <- function(conns=conns, # A data frame of co-accessibility scores, also the output of Step.3
                    GPTab=GPTabFilt,  # A list of enhancer cluster, the output of Step.2
                    cutoff=0.1 # The cutoff of co-accessibility score to determine whether the enhancer pairs are significantly co-accessible.
){
  nodes <- as.vector(GPTab$Peak)
  conn.sub <- conns %>% filter(Peak1 %in% nodes & Peak2 %in% nodes) %>% filter(coaccess >= cutoff) 
  conn.sub <- conn.sub %>% setNames(c("from","to","weight"))
  ### If the graph has no edges, return NULL
  if (nrow(conn.sub) == 0 || ncol(conn.sub) == 0) {
    return(NULL)
  }
  g <- igraph::graph_from_data_frame(d = conn.sub, vertices = data.frame(name = nodes), directed = FALSE) %>% 
    simplify(remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "mean")
  return(g)
}

# calculate several metrics for each enhancer network
NetworkComplexity <- function(net=net,
                              gene = gene,
                              GPTab=GPTabFilt)
  {
  if(class(net) == "NULL"){
    data <- data.frame(gene=gene,
                       EdgeNum=0, 
                       NetworkSize=nrow(GPTab),
                       MaxDegree=0, 
                       NetworkConnectivity=0, 
                       row.names = gene)
  }else{
    EdgeNum <- length(E(net))*2
    NetworkSize <- length(V(net))
    MaxDegree <- max(degree(net, v = V(net), mode = 'all'))
    NetworkConnectivity <- EdgeNum/NetworkSize
    data <- data.frame(gene=gene,
                       EdgeNum=EdgeNum, 
                       NetworkSize=NetworkSize,
                       MaxDegree=MaxDegree, 
                       NetworkConnectivity=NetworkConnectivity, 
                       row.names = gene)

  }
  
  return(data)
}


#' Calculate Genome-wide Conservation Score for Position-Specific Regions
#'
#' This function calculates the phastCons conservation scores for a given set of genomic regions.
#' It checks if the specified reference dataset is installed, installs it if necessary, and then computes the scores.
#'
#' @param GR A GRanges object representing the genomic regions of interest.
#' @param ref A string specifying the name of the phastCons dataset to use. Default is "phastCons30way.UCSC.hg38".
#' @return A GRanges object with an added metadata column `phastCons_Score` representing the conservation scores.
#' @import GenomicScores
#' @importFrom BiocManager install
#' @export
phastCons_Score <- function(GR, ref = "phastCons30way.UCSC.hg38") {
    # Check if the specified reference dataset is installed, install if not
    if (!requireNamespace(ref, quietly = TRUE)) {
        BiocManager::install(ref)
    }
    
    # Load the GenomicScores package and get the conservation scores
    phast <- GenomicScores::getGScores(ref)
    
    # Calculate the phastCons scores for the provided genomic regions
    # Rename the default score column to 'phastCons_Score'
    score <- GenomicScores::gscores(phast, GR) %>% mutate(phastCons_Score = default) %>% select(-default)
    
    return({score})
}
