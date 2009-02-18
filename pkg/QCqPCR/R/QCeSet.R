plotPairs <- function(x, y, eset, saveToFile) # plots graph between the 2 samples
{
    if (saveToFile) jpeg(filename <- paste(x, y, ".jpeg", sep = ""))
    raw <- as.data.frame(exprs(eset))
    plot(raw[, x], raw[, y], xlab = x, ylab = y, xlim = c(1, max(raw[, x], na.rm=T)), ylim = c(1, max(raw[, y], na.rm=T)))
    title(main <- c(x, "vs", y, "R^2 = ", cor(raw[, x] ,raw[, y], use <- "complete.obs")))
    abline(0, 1)
    if (saveToFile) dev.off()
}

pairsqPCR <- function(eset,saveToFile=FALSE)
{
for (i in sampleNames(eset))
    for (j in sampleNames(eset)) plotPairs(i, j, eset, saveToFile)
}

corMatrixqPCR <- function(eset)
{
matrix_maker <- vector()
raw <- as.data.frame(exprs(eset))
no_samples <- length(sampleNames(eset))
for (i in sampleNames(eset))
    for (j in sampleNames(eset)) matrix_maker <- c(matrix_maker,cor(raw[, i],raw[, j],use="complete.obs"))
correlation_matrix <- matrix(matrix_maker,no_samples,no_samples)
colnames(correlation_matrix) <- gsub(" ", "", sampleNames(eset))
row.names(correlation_matrix) <- gsub(" ", "", sampleNames(eset))
correlation_matrix[upper.tri(correlation_matrix)] <- NA
return(correlation_matrix)
}
