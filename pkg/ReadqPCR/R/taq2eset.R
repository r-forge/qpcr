taq2eset <- function(filename){
    raw <- read.delim(filename)
    levels(raw$Sample) <- gsub(" ","_",levels(raw$Sample))
    samples <- levels(raw$Sample)
    detectors <- levels(raw$Detector)

    eset <- data.frame(detectors,row.names=1)
    for (sample in samples) {
       Cts <- data.frame(raw$Detector[raw$Sample == sample], 
                         as.numeric(as.character(raw$Ct[raw$Sample == sample])), 
                         row.names=1)
       eset <- data.frame(merge(eset,Cts, by="row.names"), row.names=1)
    }
    names(eset) <- samples
    eset <- new("ExpressionSet", exprs=as.matrix(eset)) # make it an expression set

    return(eset)
}
