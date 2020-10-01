# runMany.R
#----------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#----------------------------------------------------------------------------------------------------
if(!interactive()) {
   args <- commandArgs(trailingOnly=TRUE)
   if(length(args) != 2){
      printf("usage:  Rscript runMany.R  <startGeneNumber> <endGeneNumber>")
      stop()
   }
   #args<-c(1,5)
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integer(args[2])
   }
#------------------------------------------------------------------------------
library(RUnit)
library(trenaSGM)
library(MotifDb)
library(BiocParallel)
library(futile.logger)
library(RPostgreSQL)
library(org.Hs.eg.db)
#------------------------------------------------------------------------------
stopifnot(packageVersion("trena") >= "1.5.13")
stopifnot(packageVersion("trenaSGM") >= "0.99.76")
#------------------------------------------------------------------------------
# we need a configuration file (of R commands), specified on the command line
# informing us of how to run this whole genome parallelized script
#------------------------------------------------------------------------------
configurationFile <- "configPlacenta_New.R" #Need to change to WD
stopifnot(file.exists(configurationFile))

if(!exists("configurationFileRead") || !configurationFileRead)
  source(configurationFile)
#------------------------------------------------------------------------------
if(!interactive()){
   args <- commandArgs(trailingOnly=TRUE)
   stopifnot(length(args) == 2)
   startGeneIndex <- as.integer(args[1])
   endGeneIndex <- as.integer(args[2])
   }
#-----------------------------------------------------------------------------
stopifnot(exists("trenaProject"))
stopifnot(exists("mtx"))
stopifnot(exists("goi"))
#OUTPUTDIR <- "~/Dropbox/PlacentalTRNProject/Placental_TRN_July/TrenaGSMTest/"
#LOGDIR <- "~/Dropbox/PlacentalTRNProject/Placental_TRN_July/TrenaGSMTest/"
if(!file.exists(OUTPUTDIR)) dir.create(OUTPUTDIR)
if(!file.exists(LOGDIR)) dir.create(LOGDIR)
#----------------------------------------------------------------------------------------------------
basic.build.spec <- list(title="testPlacenta_Snorkel",
                         type="footprint.database",
                         stageDirectory=OUTPUTDIR,#Config
                         genomeName="hg38",
                         matrix=mtx,#Config
                         db.host=getFootprintDatabaseHost(trenaProject),
                         db.port=getFootprintDatabasePort(trenaProject),
                         databases=desired.footprint.databases,#Config
                         annotationDbFile=dbfile(org.Hs.eg.db),
                         motifDiscovery="builtinFimo",
                         tfPool=tf.pool,#Config
                         tfMapping=c("MotifDB", "TFClass"),
                         #tfPrefilterCorrelation=tfPrefilterCorrelation,
                         tfPrefilterCorrelation=tfPrefilterCorrelation, #prefiltering
                         correlationThreshold=correlationThreshold, #after to restrict size
                         orderModelByColumn="pearson",
                         solverNames=SOLVERS,
                         solvers=SOLVERS) #Config

#----------------------------------------------------------------------------------------------------
buildModel <- function(short.spec)
{
  required.fields <- c("targetGene", "runParallel")
  missing.fields <- setdiff(required.fields, names(short.spec))

  if(length(missing.fields) > 0){
    msg <- sprintf("buildModel finds fields missing in short.spec: %s", paste(missing.fields, collapse=", "))
    stop(msg)
  }

  printf("********* [building model for %s, parallel = %s]", short.spec$targetGene, short.spec$runParallel)


  spec <- basic.build.spec
  targetGene <- short.spec$targetGene

  filename <- sprintf("%s/%s.RData", OUTPUTDIR, targetGene)
  system(sprintf("touch %s", filename))  # so an empty file exists if the model building fails

  tbl.geneLoc <- subset(tbl.geneInfo, geneSymbol==targetGene)[1,]
  chromosome <- tbl.geneLoc$chrom
  tss <- tbl.geneLoc$tss
  genomeName <- spec$genomeName

  spec <- basic.build.spec
  spec$targetGene <- targetGene
  spec$tss=tss
  spec$regions <- determineRegulatoryRegions(targetGene)

  spec$geneSymbol <- targetGene

  goodEnoughCor <- checkCor(targetGene, mtx)
  if(!goodEnoughCor)
    results <- list(model=data.frame(), regulatoryRegions=data.frame())
  else{
    builder <- FootprintDatabaseModelBuilder(genomeName, targetGene, spec, quiet=FALSE)
    results <- build(builder)
  }
  # to not save the footprints, go into results and make the footprints an empty dataframe.
  save(results, file=filename)

  return(results)

} # buildModel
#----------------------------------------------------------------------------------------------------
# small check to see if there are enough TFs that correlate to run a model


checkCor <- function(targetGene, mtx)
{
  target.gene.expression <- mtx[targetGene,]

  other.tfs <- intersect(rownames(mtx), allKnownTFs())

  # is our target.gene itself a TF? If so, eliminate it
  target.gene.among.tfs <- match(targetGene, other.tfs)
  if(!is.na(target.gene.among.tfs))
    other.tfs <- other.tfs[-target.gene.among.tfs]

  mtx.test <- mtx[other.tfs,]

  correlations <- abs(apply(mtx.test, 1, function(row) cor(target.gene.expression, row)))
  range.of.correlations <- fivenum(correlations)
  print(fivenum(correlations))
  third.quartile <- range.of.correlations[4]
  max <- range.of.correlations[5]

  return(third.quartile > 0)
}
#----------------------------------------------------------------------------------------------------
do.run <- function(genes, parallel=TRUE)
{
  short.specs <- lapply(genes, function(gene)
    list(targetGene=gene,
         regionsMode="enhancers",
         runParallel=parallel))

  names(short.specs) <- as.character(genes)

  if(parallel){
    bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=LOGDIR, threshold="INFO", workers=WORKERS)
    printf("running bplapply on %d genes, %d workers", length(genes), WORKERS)
    results <- bptry({bplapply(short.specs, buildModel, BPPARAM=bp.params)})
  } else {
    results <- lapply(short.specs, buildModel)
    names(results) <- genes
  }

  invisible(results)

} # do.run
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Test first 4 genes in model
test_fourGenes <- function(useParallel=FALSE)
{
   genes <- c("TSPAN6","TNMD","DPM1","SCYL3","C1orf112")
   x <- do.run(genes, parallel=useParallel)
dim(x$TSPAN6$model)[1]>0
dim(x$TNMD$model)[1]>0
dim(x$SCYL3$model)[1]>0
dim(x$C1orf112$model)[1]>0
} # test_fourGenes
#----------------------------------------------------------------------------------------------------
if(!interactive()){

   stopifnot(startGeneIndex < length(goi))
   stopifnot(endGeneIndex <= length(goi))
   stopifnot(startGeneIndex < endGeneIndex)

   goi.thisRun <- goi[startGeneIndex:endGeneIndex]
   printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
   x <- do.run(goi.thisRun, parallel=TRUE)
   }

if(interactive()){  # for development and debugging
   startGeneIndex <- 1
   endGeneIndex <- 10
   goi.thisRun <- goi[startGeneIndex:endGeneIndex]
   printf("running with genes %d - %d", startGeneIndex, endGeneIndex)
   x <- do.run(goi.thisRun, parallel=FALSE)
   }


#------------------------------------------------------------------------------
