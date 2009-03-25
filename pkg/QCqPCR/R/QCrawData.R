plotPairs <- function(x, y, qPCRSet, saveToFile) # plots graph between the 2 samples
{
    if (saveToFile) jpeg(filename <- paste(x, y, ".jpeg", sep = ""))
    raw <- as.data.frame(exprs(qPCRSet))
    plot(raw[, x], raw[, y], xlab = x, ylab = y, xlim = c(1, max(raw[, x], na.rm=T)), ylim = c(1, max(raw[, y], na.rm=T)))
    title(main <- c(x, "vs", y, "R^2 = ", cor(raw[, x] ,raw[, y], use <- "complete.obs")))
    abline(0, 1)
    if (saveToFile) dev.off()
}

pairsqPCR <- function(qPCRSet, saveToFile=FALSE)
{
for (i in sampleNames(qPCRSet))
    for (j in sampleNames(qPCRSet)) plotPairs(i, j, qPCRSet, saveToFile)
}

corMatrixqPCR <- function(qPCRSet)
{
matrix_maker <- vector()
raw <- as.data.frame(exprs(qPCRSet))
no_samples <- length(sampleNames(qPCRSet))
for (i in sampleNames(qPCRSet))
    for (j in sampleNames(qPCRSet)) matrix_maker <- c(matrix_maker, cor(raw[, i], raw[, j], use="complete.obs"))
correlation_matrix <- matrix(matrix_maker,no_samples, no_samples)
colnames(correlation_matrix) <- gsub(" ", "", sampleNames(qPCRSet))
row.names(correlation_matrix) <- gsub(" ", "", sampleNames(qPCRSet))
correlation_matrix[upper.tri(correlation_matrix)] <- NA
return(correlation_matrix)
}
