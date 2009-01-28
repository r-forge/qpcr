pairsqPCR <- function(eset,saveToFile=FALSE)
{
for (i in sampleNames(eset))
        for (j in sampleNames(eset)) plotPairs(i, j, eset, saveToFile)
}
