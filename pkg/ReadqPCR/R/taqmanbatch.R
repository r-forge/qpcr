setClass("qPCRSet", contains = "eSet")

## accessor methods
setMethod("exprs", signature = "qPCRSet", definition = 
    function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "qPCRSet", definition = 
    function (object, value) assayDataElementReplace(object, "exprs", value)
)

read.taqman <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    auxnames <- unlist(list(...))
    filenames <- c(filenames, auxnames)
    checkValidTaqmanFilenames(filenames)
    pdata <- pData(phenoData) # number of filesi
    taqInfo <- .read.TaqBatch(filenames, verbose) # need to make this work for tech reps and multiple files
    exprs <- taqInfo$exprs
    original.order <- taqInfo$origOrder
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs), extended = TRUE)
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    return(new("qPCRSet", exprs = exprs, phenoData = phenoData, well.order = original.order))
}

.read.TaqBatch <- function(filenames, verbose)
{
    for (filename in filenames) {
        if (verbose) cat("Filename[ ",filename," ]")
        raw <- read.delim(filename, skip = 12) # read in file, ignoring the gubbins on the first 12 lines
        if (! 1 %in% regexpr("Summary", as.character(raw[,1]))) stop("Problems with Taqman file, Summary info not found") # this could be made safer..ie if a sample was named Summary also its hard to understand...any better ideas anyone?
        EndOfData <- grep("Summary", raw[,1])
        raw <- raw[1:EndOfData-1, ]  # get rid of from where data finishes until the end of the file
        raw$Sample = factor(raw$Sample) # clean up additional levels brought in by extra info in the raw file before chopping
        raw$Detector = factor(raw$Detector)
#        levels(raw$Sample) <- gsub(" ", "_", levels(raw$Sample)) # replace spaces with _ for sample names USE MAKENAMES
#        levels(raw$Detector) <- gsub(" ", "_", levels(raw$Detector)) # replace spaces with _ for detectors USE MAKENAME
        levels(raw$Detector) <- make.names(levels(raw$Sample))
        levels(raw$Detector <- make.names(levels(raw$Detector))
        samples <- levels(raw$Sample)
        detectors <- levels(raw$Detector)
        original.order <- list() # initialise the list
        exprs <- data.frame(detectors, row.names=1) # start the exprs data frame

	raw$Ct[as.character(raw$Ct) %in% "Undetermined"] <- NA

        for (sample in samples) { # for each sample
            if (verbose) cat("Now reading for sample:", sample, "\n")
            # work out if there are technical replicates
            total.detectors <- length(raw$Detector[raw$Sample == sample])
            individual.detectors <- length(levels(raw$Detector[raw$Sample == sample]))
            tec.reps <- total.detectors/individual.detectors
            if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
                warning.text = paste("Corrupt taqman file: total number of readings for sample ", 
                       sample, " not a multiple of number of individual number of detectors")
                stop(warning.text)
            }
            if (tech.reps > 1) { # Currently can't cope with technical replicates
                warning.text = "More than 1 technical replicate detected"
                stop(warning.text)
            }
#            raw$Ct[as.character(raw$Ct)[raw$Sample == sample] %in% "Undetermined"] = NA # change Undertermined values to NA to stop warning messages appearing when Ct vector coerced into numeric vector
            original.order = c(original.order,list(cbind(as.character(raw$Detector[raw$Sample == sample]),
                as.character(raw$Ct[raw$Sample == sample])))) # This bit to add the information about pipetting and order

            Cts <- data.frame(raw$Detector[raw$Sample == sample], # put Cts values in a matrix
                         as.numeric(as.character(raw$Ct[raw$Sample == sample])),
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
    taqInfo$origOrder <- original.order
    return(taqInfo)
}
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
