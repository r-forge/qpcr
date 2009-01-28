plotPairs <- function(x, y, eset, saveToFile) # plots graph between the 2 samples
{
if (saveToFile) jpeg(filename <- paste(x, y, ".jpeg", sep = ""))
raw <- as.data.frame(exprs(eset))
plot(raw[, x], raw[, y], ylim <- c(1, 40), xlim <- c(1, 40), xlab <- x, ylab = y)
title(main <- c(x, "vs", y, "R^2 = ", cor(raw[, x] ,raw[, y], use <- "complete.obs")))
abline(0, 1)
if (saveToFile) dev.off()
}

