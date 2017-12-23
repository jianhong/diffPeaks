#' Call peaks
#' @description Calling peaks from bam files
#' @param bamfile bam file name
#' @param index index file name
#' @param txdb TxDb object
#' @param genome BSgenome object
#' @param upstream peak scan start position from upstream of gene
#' @param downstream peak scan end position till downstream of gene
#' @param ideaPeakWidth ideally peak width, eg, ATAC-seq: 200bp
#' @param FDRfilter cut off value of FDR.
#' @param direction over or less. In most case, default over is OK. 
#'        In some cases, such as cohesion, less may be what you want to try.
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicFeatures
#' @import BSgenome
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import ChIPpeakAnno
#' @import GenomeInfoDb
#' @importFrom stats pnbinom
#' @export
#' @return a GRanage object indicates the peaks. 
#' @examples 
#' path <- system.file("extdata", package = "diffPeaks", mustWork = TRUE)
#' bamfiles <- dir(path, "bam$", full.names = TRUE)
#' bamfile <- bamfiles[1]
#' index <- bamfile
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome <- Hsapiens
#' upstream <- downstream <- 50000L
#' ideaPeakWidth <- 200
#' direction <- "over"
#' x <- callPeak(bamfile, txdb=txdb, genome=genome)
callPeak <- function(bamfile, index=bamfile,
                     txdb, genome, 
                     upstream=50000, downstream=50000,
                     ideaPeakWidth=200, FDRfilter=0.05, 
                     direction=c("over", "less")){
  stopifnot(is(txdb, "TxDb"))
  stopifnot(is(genome, "BSgenome"))
  stopifnot(is(upstream, "numeric"))
  stopifnot(is(downstream, "numeric"))
  stopifnot(is(ideaPeakWidth, "numeric"))
  upstream <- upstream[1]
  downstream <- downstream[1]
  ideaPeakWidth <- ideaPeakWidth[1]
  direction <- match.arg(direction)
  suppressWarnings({
    genes <- genes(txdb)
    genes.ups <- promoters(genes, upstream = upstream, downstream = 0)
    genes.dws <- genes
    start(genes.dws[strand(genes)=="+"]) <- end(genes.dws[strand(genes)=="+"])
    end(genes.dws[strand(genes)=="-"]) <- start(genes.dws[strand(genes)=="-"])
    genes.dws <- promoters(genes.dws, upstream = 0, downstream = downstream)
    region <- reduce(c(genes, genes.ups, genes.dws))
    region <- trim(region)
  })
  if(length(region)<1){
    stop("Can not get peak calling region from txdb")
  }
  seqinfo <- scanBamHeader(bamfile, index=index)[[1]]$targets
  region <- region[seqnames(region) %in% names(seqinfo)]
  seqlevels(region) <- intersect(seqlevels(region), names(seqinfo))
  ## check SE or PE, if SE, estimate the fragment size and shift the reads
  pe <- testPairedEndBam(bamfile, index)
  if(!pe){
    ## estimate the fragment size
    fragmentSize <- estFragmentLength(bamfile, index, plot = FALSE)
    halfD <- floor(fragmentSize/2)
    gal <- readGAlignments(bamfile, index = index, 
                           param=ScanBamParam(scanBamFlag(isPaired = FALSE, 
                                                          isUnmappedQuery = FALSE,
                                                          isSecondaryAlignment = FALSE,
                                                          isNotPassingQualityControls = FALSE),
                                              what=scanBamWhat(),
                                              which = region))
    gal <- as(gal, "GRanges")
    mcols(gal) <- NULL
    gal <- split(gal, strand(gal))
    gal <- gal[c("+", "-")]
    gal[["+"]] <- shift(gal[["+"]], shift = halfD)
    gal[["-"]] <- shift(gal[["-"]], shift = -halfD)
    gal <- unlist(gal)
  }else{
    gal <- readGAlignmentPairs(bamfile, index = index,  
                               param=ScanBamParam(scanBamFlag(isProperPair = TRUE, 
                                                              isUnmappedQuery = FALSE,
                                                              isSecondaryAlignment = FALSE,
                                                              isNotPassingQualityControls = FALSE),
                                                  what=scanBamWhat(),
                                                  which = region))
    gal <- as(gal, "GRanges")
    mcols(gal) <- NULL
  }
  gal <- sort(gal)
  gal <- promoters(gal, upstream = 0, downstream = 1) ## only use 5 ends for peak calling
  cvg <- coverage(gal)
  cvg <- cvg[sapply(cvg, sum)>0]
  ## stplit gerions to ideaPeakWidth
  region <- region[seqnames(region) %in% names(cvg)]
  region <- tile(region, width = ideaPeakWidth)
  region <- unlist(region)
  region <- split(region, seqnames(region))
  region <- region[names(cvg)]
  counts <- viewSums(Views(cvg, region))
  region <- unlist(region, use.names = FALSE)
  region$counts <- unlist(counts, use.names = FALSE)
  region <- region[region$counts>0]
  if(length(region)==0){
    return(NULL)
  }
  ## Negative Binomial Distribution
  ## P(X=r) = C(n-1, r-1) p^r (1-p)^(n-r)
  mu <- mean(region$counts)
  var <- var(region$counts)
  p <- 1 - mu / var
  if(p < 0){
    p <- 1
  }
  r <- region$counts
  n <- width(region)
  region$pval <- 
    pnbinom(r-1, n-1, 
            prob = p, 
            lower.tail = FALSE) #choose(n-1, r-1) * p^r * (1-p)^(n-r)
  region$FDR <- p.adjust(region$pval, method = "BH")
  region <- region[region$FDR < FDRfilter]
  region.rd <- reduce(region, min.gapwidth = ideaPeakWidth, with.revmap=TRUE)
  region.revmap <- region[unlist(region.rd$revmap)]
  region.revmap$newID <- rep(seq_along(region.rd), lengths(region.rd$revmap))
  mc <- as.data.frame(mcols(region.revmap))
  fa <- formatC(mc$newID, width = nchar(max(mc$newID)), flag = "0")
  region.rd$counts <- as.numeric(rowsum(mc$counts, fa))
  region.rd$pvalue <- sapply(split(mc$pval, fa), min)
  region.rd$FDR <- sapply(split(mc$FDR, fa), min)
  region.rd$revmap <- NULL
  region.rd[ifelse(direction=="over", region.rd$counts>mu, region.rd$counts<mu)]
}
