# creates aggregate of all models in tbl.models, file is tbl.models.all.RData
# adding hugo geneSymbols for each tf, and for the target gene
#-------------------------------------------------------------------------------------------------------
library(trenaSGM)
#-------------------------------------------------------------------------------------------------------
target.dir <- "/proj/price4/apaquett/TRENA_GSMs/GSM_ForHoldout_11162019/"

stopifnot(file.exists(target.dir))
files <- grep(".RData$", list.files(target.dir), value=TRUE)
length(files)
tbls.all <- list()

for(file in files){
   full.path <- file.path(target.dir, file)
   if(file.size(file.path(target.dir, file)) == 0) next;
   target.gene <- sub(".RData", "", file)
   x <- tryCatch({
        get(load(full.path))
        }, error=function(e) {
             return(list(model=data.frame(), regulatoryRegions=data.frame()))
             }
        )
   tbl.model <- x$model
   tbl.reg <- x$regulatoryRegions
   printf("model rows: %d", nrow(tbl.model))
   if(nrow(tbl.model) == 0) next
   tbl.model$targetGene <- target.gene
   #trimmed <- trimModel(tbl.model, tbl.reg, votesNeeded=1)
   tbl <- tbl.model
   tbl$rank <- seq_len(nrow(tbl))
   printf("model for %s: %d rows, %d cols", target.gene, nrow(tbl), ncol(tbl))
   tbls.all[[target.gene]] <- tbl
   } # for file

tbl.models <- do.call(rbind, tbls.all)
dim(tbl.models)
rownames(tbl.models) <- NULL
save(tbl.models, file=file.path(target.dir, "GSM_Placenta11152019_ForHoldout.RData"))

