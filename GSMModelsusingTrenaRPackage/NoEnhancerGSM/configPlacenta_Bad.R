library(RUnit)
library(TrenaProjectPlacenta)

printf("--- reading config.R")

stopifnot(packageVersion("TrenaProject") >= "0.99.37")
#stopifnot(packageVersion("TrenaProjectBrainCell") >= "0.99.04")

trenaProject <- TrenaProjectPlacenta()

# I added some minor parcing of the expression matrix to get rid of genes with no expression.
#matrix.name <- "Micro_TYROBP"
#matrix.name <- "RNAseqData"
matrix.name<-"TrainingDataForServer"

stopifnot(matrix.name %in% getExpressionMatrixNames(trenaProject))
mtx <- getExpressionMatrix(trenaProject, matrix.name)
#mtx.df <- as.data.frame(mtx)
#mtx.df$mean <- rowMeans(mtx.df)
#mtx.df <- mtx.df[(mtx.df$mean > 0.05),]
#mtx.df <- mtx.df[!(mtx.df$mean == 0),]
#mtx.df$mean <- NULL
#mtx <- as.matrix(mtx.df)
dim(mtx)== c(14718,307)

tbl.geneHancer <- get(load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData")))

# there might be a bettter solution, but for now, I'm going to make sure there's only one gene symbol per gene here
# the reason is that the determineRegulatoryRegions function doesn't currently handle things well when there's more than one entry for tbl.geneInfo
tbl.geneInfo <- get(load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData")))
tbl.geneInfo <- tbl.geneInfo[!duplicated(tbl.geneInfo$geneSymbol),]

#OUTPUTDIR <- "~/Dropbox/PlacentalTRNProject/Placental_TRN_July/TrenaGSMTest/"
OUTPUTDIR <- "/proj/price4/apaquett/TRENA_GSMs/GSM4BAD_Placenta"
if(!file.exists(OUTPUTDIR))
   dir.create(OUTPUTDIR)

WORKERS <- 15 # changed from 20 for now
SOLVERS <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman", "sqrtLasso")

LOGDIR <- file.path(OUTPUTDIR, "logs")

trenaProject <- TrenaProjectPlacenta()


desired.footprint.databases <- getFootprintDatabaseNames(trenaProject)[7]
print(desired.footprint.databases)

### CHANGE FOR TRENA GSM CODE #########


#------------------------------------------------------------------------------------------------------------------------
# we want one gene for which we have genehancer regions, and one without.
# we also want these genes to exhibit high variance across samples in the expression matrix
# these provide good (though not exhaustive) tests of building a gsm:
#   - different method for
pickGuineaPigGenes <- function(mtx)
{
   load(system.file(package="TrenaProject", "extdata", "genomeAnnotation", "geneHancer.v4.7.allGenes.RData"))
   tbl.geneHancer <- tbl.enhancers
   dim(tbl.geneHancer)

   load(system.file(package="TrenaProject", "extdata", "geneInfoTable_hg38.RData"))
   geneSymbols.in.mtx <- rownames(mtx)
   printf("geneSymbols all chromosomes: %d", length(geneSymbols.in.mtx))

   geneSymbols.with.enhancers <- intersect(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol))
   printf("genes with geneHancer annotation: %d/%d", length(geneSymbols.with.enhancers), length(geneSymbols.in.mtx))
   geneSymbols.without.enhancers <- setdiff(geneSymbols.in.mtx, unique(tbl.geneHancer$geneSymbol))  # 3174

     # choose one maximally variable gene from each set (with and w/o enhancers)
   x <- as.numeric(lapply(geneSymbols.without.enhancers,  function(gene) sd(mtx[gene,])))
  # big.enough.sd <- fivenum(x)[1]
   #genes.with.big.enough.sd <- rownames(mtx)[which(x > big.enough.sd)]
  # candidates <- intersect(genes.with.big.enough.sd, geneSymbols.without.enhancers)   # 63
      # now make sure they are protein-coding genes
   also.protein.coding <- intersect(geneSymbols.without.enhancers, subset(tbl.geneInfo, type=="protein_coding")$geneSymbol) # 8
      # choose INKA2
  # guineaPigGenes.without.enhancers <- "INKA2"

   x <- lapply(also.protein.coding,  function(gene) sd(mtx[gene,]))
   biggest.sd <- which(as.numeric(x) == max(as.numeric(x)))
   guineaPigGenes.without.enhancers <- also.protein.coding[biggest.sd]


   x <- lapply(geneSymbols.with.enhancers,  function(gene) sd(mtx[gene,]))
   biggest.sd <- which(as.numeric(x) == max(as.numeric(x)))
   guineaPigGene.with.enhancers <- geneSymbols.with.enhancers[biggest.sd]        # GSTM1

   return(list(gene.with.enhancers=guineaPigGene.with.enhancers,
               gene.without.enhancers=guineaPigGenes.without.enhancers))

} # pickGuineaPigGenes
#------------------------------------------------------------------------------------------------------------------------
test_pickGuineaPigGenes <- function()
{
   printf("--- test_pickGuineaPigGenes")
   gpg <- pickGuineaPigGenes(mtx)
   gene.no <- gpg$gene.without.enhancers
   gene.yes <- gpg$gene.with.enhancers

      # check to see that with/without genehancer info is accurate
   checkTrue(gene.yes %in% tbl.geneHancer$geneSymbol)
   checkTrue(!gene.no %in% tbl.geneHancer$geneSymbol)

   checkTrue(sd(mtx[gene.yes,]) > 2)
   checkTrue(sd(mtx[gene.no,]) > 0.5)

} # test_pickGuineaPigGenes
#------------------------------------------------------------------------------------------------------------------------
# METHOD 2: NO GENEHANCER REGIONS
# if there are no enhancers entry: use tss +/- 5kb
# if there are enhancers, use those enhancers and +/- 2kb
#
determineRegulatoryRegions <- function(gene)
{
  tbl.concise <- tbl.geneInfo[which(tbl.geneInfo$geneSymbol == gene), c("chrom", "tss")]
  # no need to figure strand since we go 2500bp in both directions
genes<-rownames(mtx)

  if(gene %in% tbl.geneHancer$geneSymbol){
    shoulder <- 5000
    tbl.regions <- data.frame(chrom=tbl.concise$chrom,
                              start=tbl.concise$tss - shoulder,
                              end=tbl.concise$tss + shoulder,
                              stringsAsFactors=FALSE)
  }

  #setTargetGene(trenaProject, gene)
  #tbl.enhancers <- getEnhancers(trenaProject)[, c("chrom", "start", "end")]
  #shoulder <- 2000
  # tbl.promoter <- data.frame(chrom=tbl.concise$chrom,
  #                           start=tbl.concise$tss - shoulder,
  #                           end=tbl.concise$tss + shoulder,
  #                          stringsAsFactors=FALSE)
  #tbl.regions <- rbind(tbl.promoter, tbl.enhancers)
  #}

  if(!gene %in% tbl.geneHancer$geneSymbol){
    shoulder <- 5000
    tbl.regions <- data.frame(chrom=tbl.concise$chrom,
                              start=tbl.concise$tss - shoulder,
                              end=tbl.concise$tss + shoulder,
                              stringsAsFactors=FALSE)
  }

  new.order <- order(tbl.regions$start, decreasing=FALSE)
  tbl.regions <- tbl.regions[new.order,]
  rownames(tbl.regions) <- NULL
  tbl.reduced <- as.data.frame(union(GRanges(tbl.regions), GRanges(tbl.regions)))[, c("seqnames", "start", "end")]
  colnames(tbl.reduced) <- c("chrom", "start", "end")
  tbl.reduced$chrom <- as.character(tbl.reduced$chrom)

  return(tbl.reduced)

} # determineRegulatoryRegions


### GOTT HERE #####
#------------------------------------------------------------------------------------------------------------------------
test_determineRegulatoryRegions <- function()
{
   printf("--- test_determineRegulatoryRegions")

   tbl.gh <- determineRegulatoryRegions("RPS4Y1")   # has genehancer regions
   checkTrue(nrow(tbl.gh) >= 10)
   tss.gstm1 <- subset(tbl.geneInfo, geneSymbol=="RPS4Y1")$tss
   gr.10kpromoter <- GRanges(data.frame(chr="chrY", start=tss.gstm1-2000, end=tss.gstm1+2000, stringsAsFactors=FALSE))
   gr.regions <- GRanges(tbl.gh)
   checkEquals(length(findOverlaps(gr.10kpromoter, gr.regions, type="within")), 1)

   tbl.no.gh <- determineRegulatoryRegions(guineaPigGene.without.enhancers)    # no genehancer regions
   checkEquals(dim(tbl.no.gh), c(1,3))
   with(tbl.no.gh, checkEquals(end - start, 10000))

} # test_determineRegulatoryRegions
#------------------------------------------------------------------------------------------------------------------------
# generate bash script to run genes with bplapply in chunks

max <- nrow(mtx)
starts <- seq(1, max, 20)
ends <- starts + 19
ends[length(ends)] <- max
bashScript <- file("runByChunks20.sh")
lines <- sprintf("Rscript runMany.R %d %d", starts, ends)
writeLines(lines, bashScript)
close(bashScript)
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_pickGuineaPigGenes()
   test_determineRegulatoryRegions()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
#runTests()

#test.goi <- as.character(pickGuineaPigGenes(mtx))

goi <- rownames(mtx)
configurationFileRead <- TRUE
tfPrefilterCorrelation=0
correlationThreshold=0
tf.pool <- (intersect(trenaSGM::allKnownTFs(identifierType="geneSymbol"), mcols(MotifDb)$geneSymbol))
tf.pool <- intersect(tf.pool, rownames(mtx))
printf("using %d tfs, each with a MotifDb matrix and expression in mtx", length(tf.pool))
use.geneHancer <- TRUE

