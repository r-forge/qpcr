read.taqman <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    auxnames <- unlist(list(...))
    filenames <- c(filenames, auxnames)
    checkValidTaqmanFilenames(filenames)
    pdata <- pData(phenoData) # number of filesi
    exprs <- read_TaqBatch(filenames, verbose) # need to make this work for tech reps and multiple files
    n <- length(colnames(exprs))
    if (dim(pdata)[1] != n) { # so if we don't have a row for each sample in the pData matrix
        warning("Incompatible phenoData object. Created a new one using sample name data derived from raw data.\n")
        samplenames <- sub("^/?([^/]*/)*", "", colnames(exprs), extended = TRUE)
        pdata <- data.frame(sample = 1:length(samplenames), row.names = samplenames)
        phenoData <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "arbitrary numbering",
                row.names = "sample"))
    }
    return(new("ExpressionSet", exprs = exprs, phenoData = phenoData))
}

read_TaqBatch <- function(filenames, verbose)
{
    for (filename in filenames) {
	if (verbose) cat("Filename[ ",filename," ]")
        raw <- read.delim(filename)
        levels(raw$Sample) <- gsub(" ", "_", levels(raw$Sample)) # replace spaces with _ for sample names
        levels(raw$Detector) <- gsub(" ", "_", levels(raw$Detector)) # replace spaces with _ for detectors

        samples <- levels(raw$Sample)
        detectors <- levels(raw$Detector)
        exprs <- data.frame(detectors, row.names=1)
        for (sample in samples) { # for each sample
            if (verbose) cat("Now reading for sample:", sample, "\n")
            # work out if there are technical replicates
            total_detectors <- length(raw$Detector[raw$Sample == sample])
            individual_detectors <- length(levels(raw$Detector[raw$Sample == sample]))
            tech_reps <- total_detectors/individual_detectors
            if ((tech_reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
                warning_text = paste("Corrupt taqman file: total number of readings for sample ", 
                       sample, " not a multiple of number of individual number of detectors")
                stop(warning_text)
            }
            if (tech_reps > 1) { # Currently can't cope with technical replicates
                warning_text = "More than 1 technical replicate detected"
                stop(warning_text)
            }
            Cts <- data.frame(raw$Detector[raw$Sample == sample], # put Cts values in a matrix
                         as.numeric(as.character(raw$Ct[raw$Sample == sample])),
                         row.names=1)
            exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
            if (verbose) cat("sample ", sample, "read\n")
        }
        names(exprs) <- samples
        exprs <- as.matrix(exprs)
    }
    return(exprs)
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
