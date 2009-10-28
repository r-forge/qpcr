normaliseByHKG <- function(qPCRSet, hkg, verbose = FALSE){ # takes qPCRSet and HKG
  hkg <- make.names(hkg)
  if(hkg %in% featureNames(qPCRSet) == FALSE) stop (hkg," not found in file. Ensure entered housekeeping genes appear in the file")  
  hkgCts <- as.numeric(exprs(qPCRSet[hkg,]))
  if(sum(is.na(hkgCts)) > 0) warning(hkg, " May be a bad housekeeping gene to normalise with since it did not produce a reading ", sum(is.na(hkgCts)), "times out of", length(hkgCts), ".")
#hkg <- gsub("-.+$","",hkg) # regexp to remove any rubbish from end of control gene spec
  exprs(qPCRSet) = t(t(exprs(qPCRSet)) - hkgCts)
  return(qPCRSet)
}