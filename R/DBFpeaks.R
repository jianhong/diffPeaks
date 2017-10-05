#' @title differential binding analysis
#' @description use DBF to test differential peaks
#' @param counts output of \link{countTable}
#' @param ... parameters could be passed to \link[DESeq2]{DESeqDataSet}
#' @return an list of \link[GenomicRanges]{GRanges}
#' @import DESeq2
#' @import IRanges
#' @import BiocGenerics
#' @import SummarizedExperiment
#' @importFrom stats p.adjust pgamma
#' @importFrom ecodist MRM
#' @export
#' @author Jianhong Ou
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$")
#' peaks <- dir(path, "bed$")
#' p <- mergePeaks(file.path(path, peaks))
#' colData <- DataFrame(samples=bamfiles, condition=sub(".rep..bam", "", bamfiles))
#' cnt <- countTable(p, file.path(path, bamfiles), colData)
#' DBFpeaks(cnt, design= ~condition)
#' 
DBFpeaks <- function(counts, ...){
  if(length(counts$signature)!=1){
    stop("counts must be output of countTable!")
  }
  if(counts$signature!="countTable"){
    stop("counts must be output of countTable!")
  }
  if(any(names(counts)!=c("feature", "tile.feature", "signature"))){
    stop("counts must be output of countTable!")
  }
  
  ## check the length of each peak, 
  ## case 2: use the DESeq2 results, combine two results
  ## case 3: use the DESeq2 results, combine three results
  ## case 4~: use DBF.test
  colData <- colData(counts$tile.feature)
  gr <- rowRanges(counts$tile.feature)
  comparison <- ""
  dds <- DESeqDataSet(counts$tile.feature, ...)
  dds <- DESeq(dds, betaPrior = FALSE, fitType = "local")
  rld <- rlog(dds, blind = FALSE)
  tile.cnt <- assay(rld)
  tile.cnt <- as.data.frame(tile.cnt)
  if(length(resultsNames(dds))==2 && resultsNames(dds)[1]=="Intercept"){
    res <- results(dds)
    resLFC <- lfcShrink(dds, coef = 2, res = res)
    gr1 <- gr
    mcols(gr1) <- cbind(mcols(gr), resLFC)
    comparison <- strsplit(resultsNames(dds)[2], "_")[[1]][-3]
  }else{
    stop("Only simple comparison is supported.")
  }
  stopifnot(length(comparison)==3)
  groups <- colData[, comparison[1]]
  gpA <- comparison[2]
  gpB <- comparison[3]
  if(sum(groups %in% gpA)!=sum(groups %in% gpB)){
    stop("unbalanced groups are not supported.")
  }
  tile.cnt.s <- split(tile.cnt, as.character(gr$X))
  tile.cnt.dbf <- lapply(tile.cnt.s, function(.ele){
    data <- rbind(as.matrix(unname(.ele[, which(groups %in% gpA)])), 
                  as.matrix(unname(.ele[, which(groups %in% gpB)])))
    rownames(data) <- NULL
    maxDist <- rowSums(data)
    maxDist <- maxDist[seq.int(nrow(.ele))] - 
      maxDist[nrow(.ele) + (seq.int(nrow(.ele)))]
    maxDist <- maxDist[which.max(abs(maxDist))]
    baseMean <- 2^mean(as.matrix(.ele))
    group.labels <- rep(c(gpA, gpB), each=nrow(.ele))
    stats <- tryCatch(DBF.test(as.matrix(ecodist::distance(data, "mahalanobis")), 
                               group.labels, nrow(data)),
                      error=function(e){return(c(NA, NA))})
    c(baseMean=baseMean, log2FoldChange=maxDist, stats)
  })
  tile.cnt.n <- as.numeric(names(tile.cnt.dbf))
  tile.cnt.dbf <- do.call(rbind, tile.cnt.dbf)
  dbf.res <- data.frame(X=tile.cnt.n, tile.cnt.dbf)
  dbf.res <- dbf.res[match(seq_along(rowRanges(counts$feature)), dbf.res$X), ]
  dbf.res$padj <- p.adjust(dbf.res$dbf.p.value, method="BH")
  ## case 2
  ## case 3
  gr1 <- gr1[order(gr1$pvalue, decreasing = FALSE)]
  gr1 <- gr1[!duplicated(gr1$X)]
  gr1 <- gr1[order(gr1$X, decreasing = FALSE)]
  gr1$range <- paste(start(gr1), end(gr1), sep="-")
  ranges(gr1) <- ranges(counts$feature)
  gr1$type <- "subRange"
  ## case 4~
  nNA <- !is.na(dbf.res$dbf.p.value)
  gr1[nNA]$type <- "curveComparison"
  gr1[nNA]$baseMean <- dbf.res[nNA, "baseMean"]
  gr1[nNA]$log2FoldChange <- dbf.res[nNA, "log2FoldChange"]
  gr1[nNA]$stat <- dbf.res[nNA, "dbf.statistic"]
  gr1[nNA]$pvalue <- dbf.res[nNA, "dbf.p.value"]
  gr1[nNA]$padj <- dbf.res[nNA, "padj"] ## adjust p value too high.
  gr2 <- split(gr, as.character(gr$X))
  gr2 <- unlist(GRangesList(lapply(gr2, range)))
  gr2 <- gr2[order(as.numeric(names(gr2)))]
  stopifnot(identical(as.integer(names(gr2)), seq_along(gr1)))
  gr1[nNA]$range <- paste(start(gr2), end(gr2), sep="-")[nNA]
  gr1$X <- NULL
  gr1
}