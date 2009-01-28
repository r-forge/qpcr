taq2eset <- function(filename){
    raw <- read.delim(filename)
    levels(raw$Sample) <- gsub(" ", "_", levels(raw$Sample)) # replace spaces with _ for sample names
    levels(raw$Detector) <- gsub(" ", "_", levels(raw$Detector)) # replace spaces with _ for detectors

    samples <- levels(raw$Sample)
    detectors <- levels(raw$Detector)

    eset <- data.frame(detectors, row.names=1)
    for (sample in samples) {
       # work out if there are technical replicates
       total_detectors <- length(raw$Detector[raw$Sample == sample])
       individual_detectors <- length(levels(raw$Detector[raw$Sample == sample]))
       tech_reps <- total_detectors/individual_detectors

       if ((tech_reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
          warning_text = paste("Corrupt taqman file: total number of readings for sample ",sample," not a multiple of number of individual detectors")
          warning(warning_text)
          return(0)
       }
       if (tech_reps > 1) { # Currently can't cope with technical replicates
          warning_text = "More than 1 technical replicate detected"
          warning(warning_text)
          return(0)
       }
       Cts <- data.frame(raw$Detector[raw$Sample == sample],
                         as.numeric(as.character(raw$Ct[raw$Sample == sample])),
                         row.names=1)
       eset <- data.frame(merge(eset, Cts, by="row.names"), row.names=1)
    }
    names(eset) <- samples
    eset <- new("ExpressionSet", exprs=as.matrix(eset)) # make it an expression set
    return(eset)
}
