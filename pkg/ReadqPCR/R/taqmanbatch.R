setClass("qPCRSet", contains = "eSet")

## accessor methods
setMethod("exprs", signature = "qPCRSet", definition = 
    function (object) assayDataElement(object, "exprs")
)

setReplaceMethod("exprs", signature = "qPCRSet", definition = 
    function (object, value) assayDataElementReplace(object, "exprs", value)
)

setGeneric("exprs.well.order",
    function(object)
    standardGeneric("exprs.well.order")
)

setGeneric("exprs.well.order<-",
    function(object, ..., value)
    standardGeneric("exprs.well.order<-")
)

setMethod("exprs.well.order", signature = "qPCRSet", definition = 
    function (object) assayDataElement(object, "exprs.well.order")
)

setReplaceMethod("exprs.well.order", signature = "qPCRSet", definition = 
    function (object, value) assayDataElementReplace(object, "exprs.well.order", value)
)

read.taqman <- function(..., filenames = character(0), phenoData = new("AnnotatedDataFrame"), notes = "", verbose = FALSE)
{
    auxnames <- unlist(list(...))
    filenames <- c(filenames, auxnames)
    checkValidTaqmanFilenames(filenames)
    pdata <- pData(phenoData) # number of files
    taqInfo <- .read.TaqBatch(filenames, verbose) # need to make this work for tech reps and multiple files
    exprs <- taqInfo$exprs
    well.order <- taqInfo$well.order
    exprs.well.order <- assayDataNew("environment", exprs.well.order = well.order)
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
    return(new("qPCRSet", exprs = exprs, phenoData = phenoData, exprs.well.order = well.order))
}

.read.TaqBatch <- function(filenames, verbose)
{
cat("ERE\n")
    totalPlateIds <- vector()
    plate.offset <- 0 # this is used to name plates when combining several files 
    fileNameCount <- 1
    for (filename in filenames) {
        if (verbose) cat("Filename[ ",filename," ]")
        raw.data <- read.delim(filename, skip = 12) # read in file, ignoring the gubbins on the first 12 lines
        if (! 1 %in% regexpr("Summary", as.character(raw.data[,1]))) stop("Problems with Taqman file, Summary info not found") #
        EndOfData <- grep("Summary", raw.data[,1])
        raw.data <- raw.data[1:EndOfData-1, ] # get rid of from where data finishes until the end of the file
        raw.data$Sample = factor(raw.data$Sample) # clean up additional levels brought in by extra info in the raw.data file before chopping
        raw.data$Detector = factor(raw.data$Detector) # ditto
        levels(raw.data$Sample) <- make.names(levels(raw.data$Sample)) # make sure names are legal
        levels(raw.data$Detector) <- make.names(levels(raw.data$Detector))
        samples <- levels(raw.data$Sample) # individual sample names
        detectors <- levels(raw.data$Detector) # individual detector names
######################################################################################################
# Now we have our original info on samples and detectors we can work out if there is more than one file, whether to combine horizontally or vertically 
# Also Check if we have duplicate plate IDs

        if(fileNameCount > 1) { # do some checking
          if(TRUE %in% (samples == colnames(totalExprs))) stop("Can't combine files, > 1 sample labels are the same between samples")
          else if(! FALSE %in% (sort(detectors) == sort(rownames(totalExprs)))) cat("we are combining files with the same detector names\n")
          else stop("Problem combining files on detector names. Make sure detector names match for all files\n")

          if(TRUE %in% (raw.data$PlateID %in% totalPlateIds)) stop ("Can't proceed, duplicate plate Ids in different files. All plate IDs should be unique")
        }
################################################################################################
####################################################
# Work out if we have technical replicates
####################################################
#        no.of.detectors <- length(detectors)
#        total.detectors <- length(raw.data$Detector)
#        tech.reps <- total.detectors/total.detectors
#cat("aa", diff.detectors, "bb", max.detectors, "\n")
#        if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
#            warning.text = paste("Corrupt taqman file: total number of readings for sample ", 
#               sample, " not a multiple of number of individual number of detectors")
#            stop(warning.text)
#        }
#        #staticDetector <- raw.data$Detector
#        if (tech.reps > 1) { # Currently can't cope with technical replicates
#            warning.text = "More than 1 technical replicate detected\n"
#            cat(warning.text, "aa", diff.detectors, "bb", max.detectors)

#            for(techDetect in  unique(raw.data$Detector)) {
#                techDLength <- sum(techDetect %in% raw.data$Detector)
#                cat("LENGTH IS: ", techDLength, "\n")
#                aaaa <- paste(raw.data$Detector[raw.data$Detector %in% techDetect], 1:tech.reps, sep="_ZZ_")
#                cat("detectorNewNames are :", aaaa, "\n")
#                raw.data$Detector[raw.data$Detector %in% techDetect] <- aaaa
 #           }
 #       }
        original.order <- list() # initialise the list
        
        well.order <- data.frame(detectors, row.names=1)
	raw.data$Ct[as.character(raw.data$Ct) %in% "Undetermined"] <- NA
############################################################
# Add Plate ID information IF there were none to being with
        if ("" %in% as.character(raw.data$PlateID)) {
          wells.per.plate <- max(as.numeric(as.character(raw.data$Well)))
          number.of.plates <- length(as.character(raw.data$Well)) / wells.per.plate
          well.names <- vector(length = length(as.character(raw.data$Well)))
          plate.well <- 1
          for (plate.number in 1:number.of.plates) {
            for (well.number in 1:wells.per.plate) {
               well.names[plate.well] <- paste(plate.number + plate.offset, well.number, sep= "-")
               plate.well <- plate.well + 1
            }
          }
          plate.offset <- plate.offset + number.of.plates # add the number of plates to offset the plate number
          raw.data$PlateID <- well.names
        }
        else { ## otherwise we just paste the well number onto the Plate ID value
          totalPlateIds <- raw.data$PlateID # put them in a variable for checking for duplication
          raw.data$PlateID <- paste(raw.data$PlateID, as.character(raw.data$Well), sep= "-")
        }
        firstTimeFlag <- TRUE
################################################################
        for (sample in samples) { # for each sample
            if (verbose) cat("Now reading for sample:", sample, "\n")

            total.detectors <- length(raw.data$Detector[raw.data$Sample == sample])
            cat("total.detectors:", total.detectors, "\n")
            individual.detectors <- length(unique(raw.data$Detector[raw.data$Sample == sample]))
#if(individual.detectors * length(samples) != individual.detectors) {stop("HERE IT GOES")
            cat("individual.detectors:", individual.detectors, "\n")
            tech.reps <- total.detectors/individual.detectors
cat("Jus before\n")
            if ((tech.reps %% 1) != 0) { # if total number of replicates not a multiple of number of individual detectors
                warning.text = paste("Corrupt taqman file: total number of readings for sample ", 
                       sample, " not a multiple of number of individual number of detectors")
                stop(warning.text)
            }
            if (tech.reps > 1) { # Currently can't cope with technical replicates
                warning.text = "More than 1 technical replicate detected\n"
                cat(warning.text)
                #if(raw.data$Detector
#                cat(raw.data$Detector)
#                [raw.data$Sample == sample]
raw.data$Detector <- as.character(raw.data$Detector) # coerce to stop funny behaviour
staticDetector <- raw.data$Detector[raw.data$Sample == sample]
#cat(staticDetector, "FOR DA FIRST\n")
                 for(techDetect in unique(raw.data$Detector[raw.data$Sample == sample]))  {
#                   cat(techDetect, "for the win\n\n")
#                   cat(staticDetector, "\n\n")
                   techDLength <- sum(staticDetector %in% techDetect)
                   cat("LENGTH IS: ", techDLength, "\n")
                   aaaa <- paste(techDetect, 1:techDLength, sep="_") 
                   #aaaa <- paste(raw.data$Detector[raw.data$Detector %in% techDetect], 1:techDLength, sep="")
                   cat("detectorNewNames are :", aaaa, "\n")
                   raw.data$Detector[raw.data$Sample == sample][raw.data$Detector[raw.data$Sample == sample] %in% techDetect] <- aaaa
                   
                   #row.names(exprs) <- raw.data$Detector
#stop()
                 }
            }
            if(firstTimeFlag == TRUE) {
                exprs <- data.frame(detectors, row.names=1) # start the exprs data frame
                firstTimeFlag <- FALSE
            }

cat("Leaving this tech reps bit\n")
cat("AND NOW\n")
cat(raw.data$Detector, "raw data detectors\n")
cat("done")
raw.data$Detector <- as.factor(raw.data$Detector)
cat("now lets blow this up\n\n")
#            }
            #cat(raw.data$Detector)
            original.order = c(original.order,list(cbind(as.character(raw.data$Detector[raw.data$Sample == sample]),
                as.character(raw.data$Ct[raw.data$Sample == sample])))) # This bit to add the information about pipetting and order

            well.info <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
                         raw.data$PlateID[raw.data$Sample == sample],
                         row.names=1)
            Cts <- data.frame(raw.data$Detector[raw.data$Sample == sample], # put Cts values in a matrix
                         as.numeric(as.character(raw.data$Ct[raw.data$Sample == sample])),
                         row.names=1)
            print(Cts)
            exprs <- data.frame(merge(exprs, Cts, by="row.names"), row.names=1)
            print(exprs)
################## if we have all the same detectors per sample:
            if (! FALSE %in% (row.names(exprs) == row.names(Cts))) stop("BYE")

            well.order <- data.frame(merge(well.order, well.info, by="row.names"), row.names=1)
            if (verbose) cat("sample ", sample, "read\n")
#zzzz <- rep("\n", 100)
#cat(zzzz)
            #stop("stop it")
        }
        names(well.order) <- samples
        names(exprs) <- samples
        exprs <- as.matrix(exprs)
        well.order <- as.matrix(well.order)
        if(fileNameCount == 1) {
          totalExprs <- exprs
          totalWell.order <- well.order
        }
        else {
          totalWell.order <- cbind(totalWell.order, well.order)
          totalExprs <- cbind(totalExprs,exprs)
        }
        fileNameCount <- fileNameCount + 1
    }
    taqInfo <- list()
    taqInfo$exprs <- totalExprs
    taqInfo$well.order <- totalWell.order
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
