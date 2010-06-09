
read.qPCR <- function(filename = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    pdata <- pData(phenoData)
    checkValidqPCRFilename(filename)
    qPCRInfo <- .read.qPCR(filename, verbose) # need to make this work for tech reps and multiple files
    exprs <- qPCRInfo$exprs
    well.order <- qPCRInfo$well.order

exprs.well.order <- assayDataNew("environment", exprs.well.order = exprs)
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs))
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    if(! is.null(qPCRInfo$well.order)) {
        return(new("qPCRBatch", exprs = exprs, phenoData = phenoData, exprs.well.order = well.order))
    }
    else {

        return(new("qPCRBatch", exprs = exprs, phenoData = phenoData))
    }
}

.read.qPCR <- function(filename, verbose)
{
    noWellData <- FALSE

    raw.data <- read.table(filename, header=TRUE)
    if(is.null(raw.data$Well) || is.null(raw.data$PlateID)) {
         noWellData <- TRUE
         if (verbose) cat("No Well and/or Plate info found, skipping this part", "\n")
    }
    else {
        raw.data$PlateID <- paste(raw.data$PlateID, as.character(raw.data$Well), sep= "-")
    }
    levels(raw.data$Sample) <- make.names(levels(raw.data$Sample))
    levels(raw.data$Detector) <- make.names(levels(raw.data$Detector))
    Ct <- as.character(raw.data$Ct)
    samples <- levels(raw.data$Sample)
    detectors <- levels(raw.data$Detector)
    allDetectors <- raw.data$Detector
    firstTimeFlag <- TRUE
    for (sample in samples) { # for each sample
        if (verbose) cat("Now reading for sample:", sample, "\n")
        total.detectors <- length(allDetectors[raw.data$Sample == sample])
        individual.detectors <- length(levels(allDetectors[raw.data$Sample == sample]))
        tech.reps <- total.detectors/individual.detectors
        raw.data$Detector <- as.character(raw.data$Detector)
          if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
            warning.text = paste("File incorrect, make sure that detectors are the same for all samples")
            stop(warning.text)
        }
        if (tech.reps > 1) {
            if(verbose) cat ("More than 1 technical replicate detected\n")
              staticDetector <- raw.data$Detector[raw.data$Sample == sample]
              for(techDetect in unique(raw.data$Detector[raw.data$Sample == sample]))  {
                techDLength <- sum(staticDetector %in% techDetect)
                suffixedNames <- paste(techDetect, 1:techDLength, sep="_TechReps.")
                raw.data$Detector[raw.data$Sample == sample][raw.data$Detector[raw.data$Sample == sample] 
                  %in% techDetect] <- suffixedNames
            }
       }
       if(firstTimeFlag == TRUE) {
           exprs <- data.frame(unique(raw.data$Detector), row.names=1) # start the exprs data frame
           well.order <- data.frame(unique(raw.data$Detector), row.names=1)
           firstTimeFlag <- FALSE
       }
       raw.data$Detector <- as.factor(raw.data$Detector)
        if(noWellData == FALSE) {
            well.info <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put well info values in a matrix
              raw.data$PlateID[raw.data$Sample == sample],
                row.names=1)
        }

        Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
          as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
            row.names=1)
        exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
        if(noWellData == FALSE) {
            well.order <- data.frame(merge(well.order, well.info, by="row.names"), row.names=1)
        }
        if (verbose) cat("sample ", sample, "read\n")
    }
    qPCRInfo <- list()
    names(exprs) <- samples
    qPCRInfo$exprs <- as.matrix(exprs)
#    colnames(well.order) <- names(exprs)
    if(noWellData == FALSE) {
      colnames(well.order) <- names(exprs)
      qPCRInfo$well.order <- well.order
#      colnames(well.order) <- names(exprs)
    }
    return(qPCRInfo)
}

checkValidqPCRFilename <- function (filename)
{
    if (!is.character(filename))
        stop(strwrap(paste("file name must be specified using a character",
            "vector, not a", sQuote(typeof(filename)))), call. = FALSE)
    if (length(filename) == 0)
        stop("no file name provided")
    if (any(sapply(filename, nchar) < 1))
        stop("empty file name not allowed")
    finfo <- file.info(filename)
    whBad <- sapply(finfo[["isdir"]], function(x) !identical(FALSE,
        x))
    if (any(whBad)) {
        msg <- paste("not valid file:\n",
            paste("  ", filename[whBad], collapse = "\n"))
        stop(msg, call. = FALSE)
    }
    TRUE
}
