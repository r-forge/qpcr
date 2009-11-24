
read.qPCR <- function(filename = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    pdata <- pData(phenoData)
    checkValidqPCRFilename(filename)
    qPCRInfo <- .read.qPCR(filename, verbose) # need to make this work for tech reps and multiple files
    exprs <- qPCRInfo$exprs
    well.order <- qPCRInfo$well.order

exprs.well.order <- assayDataNew("environment", exprs.well.order = exprs)
print(str(exprs))
#    if(! is.null(qPCRInfo$well.order)) {
#     exprs.well.order <- assayDataNew("environment", exprs.well.order = well.order)
#cat("AND TO HERE\n")
#    }
#    else {
#        well.order <- matrix()
#    }
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs), extended = TRUE)
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    print(str(well.order))
    if(! is.null(qPCRInfo$well.order)) {
cat("SSS\n\n")
        return(new("qPCRSet", exprs = exprs, phenoData = phenoData, exprs.well.order = well.order))
    }
    else {
print(exprs)
print(phenoData)
cat("QQQ\n\n")

        return(new("qPCRSet", exprs = exprs, phenoData = phenoData))
    }
}

.read.qPCR <- function(filename, verbose)
{
    noWellData <- FALSE

    raw.data <- read.table(filename)
    if(is.null(raw.data$Well) || is.null(raw.data$PlateID)) {
         noWellData <- TRUE
#        if (verbose) cat("No Well and/or Plate info found, skipping this part")
cat("No Well and/or Plate info found, skipping this part\n")
    }
    else {
        raw.data$PlateID <- paste(raw.data$PlateID, as.character(raw.data$Well), sep= "-")
    }
    original.order <- list()
    levels(raw.data$Sample) <- make.names(levels(raw.data$Sample))
    levels(raw.data$Detector) <- make.names(levels(raw.data$Detector))
    Ct <- as.character(raw.data$Ct)
    samples <- levels(raw.data$Sample)
    detectors <- levels(raw.data$Detector)
#    print(raw.data)
#    cat("\n")
    well.order <- data.frame(detectors, row.names=1)
    exprs <- data.frame(detectors, row.names=1) # start the exprs data frame
    for (sample in samples) { # for each sample
        if (verbose) cat("Now reading for sample:", sample, "\n")
        total.detectors <- length(raw.data$Detector[raw.data$Sample == sample])
        individual.detectors <- length(levels(raw.data$Detector[raw.data$Sample == sample]))
        tech.reps <- total.detectors/individual.detectors
        if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
            warning.text = paste("File incorrect, make sure that detectors are the same for all samples")
            stop(warning.text)
        }
        if (tech.reps > 1) { # Currently can't cope with technical replicates

#            warning.text = "More than 1 technical replicate detected"
#            stop(warning.text)
        }
            original.order = c(original.order,list(cbind(as.character(raw.data$Detector[raw.data$Sample == sample]), as.character(raw.data$Ct[raw.data$Sample == sample])))) # This bit to add the information about pipetting and order

        if(noWellData == FALSE) {
            well.info <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
              raw.data$PlateID[raw.data$Sample == sample],
                row.names=1)
          #totalPlateIds <- raw.data$PlateID # put them in a variable for checking for duplication
        }

        Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
          as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
            row.names=1)
        exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
        if(noWellData == FALSE) {
            well.order <- data.frame(merge(well.order, well.info, by="row.names"), row.names=1)
        }
cat("AND HERE????\n")
        if (verbose) cat("sample ", sample, "read\n")
    }
    qPCRInfo <- list()
    names(exprs) <- samples
    qPCRInfo$exprs <- as.matrix(exprs)
    if(noWellData == FALSE) qPCRInfo$well.order <- well.order
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