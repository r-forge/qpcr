setClass("qPCRSet", contains = "eSet")

## accessor methods
setMethod("exprs", signature = "qPCRSet", definition = 
    function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "qPCRSet", definition = 
    function (object, value) assayDataElementReplace(object, "exprs", value)
)



#setMethod("well.order", signature = "qPCRSet", definition = 
#    function (object) assayDataElement(object, "well.order")
#)

#setReplaceMethod("well.order", signature = "qPCRSet", definition = 
#    function (object, value) assayDataElementReplace(object, "well.order", value)
#)

read.taqman <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    auxnames <- unlist(list(...))
    filenames <- c(filenames, auxnames)
    checkValidTaqmanFilenames(filenames)
    pdata <- pData(phenoData) # number of filesi
    taqInfo <- .read.TaqBatch(filenames, verbose) # need to make this work for tech reps and multiple files
    exprs <- taqInfo$exprs
    original.order <- taqInfo$origOrder
#    original.order <- assayDataNew("environment", )
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs), extended = TRUE)
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
#    return(new("qPCRSet", exprs = exprs, phenoData = phenoData, well.order = original.order))
     return(new("qPCRSet", exprs = exprs, phenoData = phenoData))
}

.read.TaqBatch <- function(filenames, verbose)
{
    for (filename in filenames) {
        if (verbose) cat("Filename[ ",filename," ]")
        raw.data <- read.delim(filename, skip = 12) # read in file, ignoring the gubbins on the first 12 lines
        if (! 1 %in% regexpr("Summary", as.character(raw.data[,1]))) stop("Problems with Taqman file, Summary info not found") # this could be made safer..ie if a sample was named Summary also its hard to understand...any better ideas anyone?
        EndOfData <- grep("Summary", raw.data[,1])
        raw.data <- raw.data[1:EndOfData-1, ]  # get rid of from where data finishes until the end of the file
        raw.data$Sample = factor(raw.data$Sample) # clean up additional levels brought in by extra info in the raw.data file before chopping
        raw.data$Detector = factor(raw.data$Detector)
#        levels(raw.data$Sample) <- gsub(" ", "_", levels(raw.data$Sample)) # replace spaces with _ for sample names USE MAKENAMES
#        levels(raw.data$Detector) <- gsub(" ", "_", levels(raw.data$Detector)) # replace spaces with _ for detectors USE MAKENAME
        levels(raw.data$Sample) <- make.names(levels(raw.data$Sample))
        levels(raw.data$Detector) <- make.names(levels(raw.data$Detector))
        samples <- levels(raw.data$Sample)
        detectors <- levels(raw.data$Detector)
        original.order <- list() # initialise the list
        exprs <- data.frame(detectors, row.names=1) # start the exprs data frame

	raw.data$Ct[as.character(raw.data$Ct) %in% "Undetermined"] <- NA
cat("here")
        wells.per.sample <- length(as.character(raw.data$Detector))/length(samples)
cat("ncol is:", length(samples), "\n")
cat("nrow is:", wells.per.sample, "\n")
        well.order <- matrix(ncol= length(samples), nrow= wells.per.sample)
        row.names(well.order) <- 1:wells.per.sample
cat("and here")

############ If we have plate IDs (ie we have no empty PlateID cells) then we name the columns these 
# if (! "" %in% as.character(raw.data$PlateID))  colnames(well.order) <- as.character(unique(raw.data$PlateID))
#k <- 1
#for (i in 1:length(samples)) {
#  for (j in 1:(wells.per.sample)) {
#    well.order[j, i] <- as.character(raw.data$Detector)[k]
#    k <- k+1
#  }
#}
print(well.order)
#well.order <- matrix(ncol= length(samples), nrow= wells.per.sample)


## Add Plate ID information
if ("" %in% as.character(raw.data$PlateID)) {
  wells.per.plate <- max(as.numeric(as.character(raw.data$Well)))
  number.of.plates <- length(as.character(raw.data$Well)) / wells.per.plate
  well.names <- vector(length = length(as.character(raw.data$Well))
  plate.well = 1
  for (plate.number in number.of.plates) {
    for (well.number in wells.per.plate) {
      well.names[plate.well] <- paste(plate.number, well.number, sep= "-")
      plate.well = plate.well + 1
    }
  }
}
#if ("" %in% as.character(raw.data$PlateID)) { # so if we find that we have empty plate IDs

}
else {

}

        for (sample in samples) { # for each sample
            if (verbose) cat("Now reading for sample:", sample, "\n")
            # work out if there are technical replicates
            total.detectors <- length(raw.data$Detector[raw.data$Sample == sample])
            individual.detectors <- length(levels(raw.data$Detector[raw.data$Sample == sample]))
            tech.reps <- total.detectors/individual.detectors
            if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
                warning.text = paste("Corrupt taqman file: total number of readings for sample ", 
                       sample, " not a multiple of number of individual number of detectors")
                stop(warning.text)
            }
            if (tech.reps > 1) { # Currently can't cope with technical replicates
                warning.text = "More than 1 technical replicate detected"
                stop(warning.text)
            }
#            raw.data$Ct[as.character(raw.data$Ct)[raw.data$Sample == sample] %in% "Undetermined"] = NA # change Undertermined values to NA to stop warning messages appearing when Ct vector coerced into numeric vector
            original.order = c(original.order,list(cbind(as.character(raw.data$Detector[raw.data$Sample == sample]),
                as.character(raw.data$Ct[raw.data$Sample == sample])))) # This bit to add the information about pipetting and order

            Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
                         as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
                         row.names=1)
            exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
            if (verbose) cat("sample ", sample, "read\n")
        }
        names(original.order) <- samples
        names(exprs) <- samples
        exprs <- as.matrix(exprs)
    }

    taqInfo <- list()

    taqInfo$exprs <- exprs
#    taqInfo$origOrder <- original.order
    return(taqInfo)
}
#calculateWellOrder <- function (f
checkValidTaqmanFilenames <- function (filenames)
{
    if (!is.character(filenames))
        stop(strwrap(paste("file names must be specified using a character",
            "vector, not a", sQuote(typeof(filenames)))), call. = FALSE)
    if (length(filenames) == 0)
        stop("no file names provided")
    if (any(sapply(filenames, nchar) < 1))
        stop("empty file names are not allowed")
    finfo <- file.info(filenames)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE,
        x))
    if (any(whBad)) {
        msg <- paste("the following are not valid files:\n",
            paste("  ", filenames[whBad], collapse = "\n"))
        stop(msg, call. = FALSE)
    }
    TRUE
}
