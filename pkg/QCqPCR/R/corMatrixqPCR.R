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
