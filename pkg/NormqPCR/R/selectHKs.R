###############################################################################
## selection of housekeepers (HKs)
###############################################################################
## x: matrix or data.frame with data from at least 3 HKs
## method: method for the computation
## minNrHKs: minimum number of HKs which should be considered
## log: logical, data on log scale?
## Symbols: character, symbols for variables/columns
## trace: locical, print information
## na.rm: remove NA values
selectHKs <- function(x, group, method = "geNorm", minNrHKs = 2, 
                      log = TRUE, Symbols, trace = TRUE, na.rm = TRUE){
    if(method == "geNorm") {
        if(class(x) == "qPCRBatch") x <- t(exprs(x))
    }
    if(method == "NormFinder") {
        if(class(x) == "qPCRBatch"){
            if (missing(group)) group <- pData(x)[,"Group"]
            x <- t(exprs(x))
        }
        else stop("'x' must be of class qPCRBatch")
    }
    if(!is.matrix(x) & !is.data.frame(x))
        stop("'x' needs to be of class matrix or data.frame")
    if(is.data.frame(x)) x <- data.matrix(x)
    n <- ncol(x)
    if(n < 3)
        stop("you need data from at least 3 variables/columns")
    if(minNrHKs >= n)
        stop("'minNrHKs' must be smaller than 'ncol(x)'")
    if(minNrHKs < 2){
        warning("'minNrHKs' < 2 => 'minNrHKs' is set to 2")
        minNrHKs <- 2
    }
    if(missing(Symbols))
        stop("'Symbols' has to be specified")
    if(length(Symbols) != n)
        stop("'Symbols' has wrong length")

    if(method == "geNorm"){
        V <- numeric(n-minNrHKs)
        names(V) <- paste(((n-1):minNrHKs), "/", (n:(minNrHKs+1)), sep = "")
        meanM <- numeric(n-minNrHKs+1)
        names(meanM) <- as.character(n:minNrHKs)
        R <- character(n)
        names(R) <- as.character(c(rep(1, minNrHKs),(minNrHKs+1):length(R)))
        for(i in n:minNrHKs){
            M <- stabMeasureM(x, log = log, na.rm = na.rm)
            names(M) <- Symbols
            ind <- which.max(M)
            meanM[n-i+1] <- mean(M)
            if(i == minNrHKs)
                R[1:minNrHKs] <- Symbols
            else
                R[i] <- Symbols[ind]

            if(i > 2){
                if(log){
                    NF.old <- rowMeans(x)
                    NF.new <- rowMeans(x[,-ind])
                    V[n-i+1] <- sd(NF.new - NF.old, na.rm = na.rm)
                }else{
                    NF.old <- apply(x, 1, geomMean, na.rm = na.rm)
                    NF.new <- apply(x[,-ind], 1, geomMean, na.rm = na.rm)
                    V[n-i+1] <- sd(log2(NF.new/NF.old), na.rm = na.rm)
                }
            }

            if(trace){
                cat("###############################################################\n")
                cat("Step ", n-i+1, ":\n")
                cat("stability values M:\n")
                print(sort(M))
                cat("average stability M:\t", meanM[n-i+1], "\n")
                if(i > 2){
                    cat("variable with lowest stability (largest M value):\t", Symbols[ind], "\n")
                    cat("Pairwise variation, (", i-1, "/", i, "):\t", V[n-i+1], "\n")
                }
            }
            x <- x[,-ind]
            Symbols <- Symbols[-ind]
        }
        return(list(ranking = R, variation = V, meanM = meanM))
    }
    if(method == "NormFinder"){

        NF <- stabMeasureRho(x, group = group, log = log, na.rm = na.rm, returnAll = TRUE)
        k <- length(NF$rho)
        R <- character(minNrHKs)
        rho <- NF$rho
        R[1] <- Symbols[which.min(rho)]
        b <- integer(minNrHKs)
        b[1] <- which.min(rho)
        rho.min <- numeric(minNrHKs)
        rho.min[1] <- rho[b[1]]

        if(trace){
            cat("###############################################################\n")
            cat("Step ", 1, ":\n")
            cat("stability values rho:\n")
            print(sort(rho))
            cat("variable with highest stability (smallest rho value):\t", Symbols[b[1]], "\n")
        }
        for(i in 2:minNrHKs){
            rho[b[i-1]] <- NA

            for (j in (c(1:k)[-b])){
                a <- c(b,j)
                a1 <- NF$d[a,]
                a2 <- colMeans(a1)*sqrt(k/(k-i))
                b1 <- NF$v[a,]
                b2 <- colMeans(b1)/i
                rho[j] <- mean(abs(a2)+sqrt(b2))
            }
            b[i] <- which.min(rho)
            R[i] <- Symbols[b[i]]
            rho.min[i] <- rho[b[i]]
            if(trace){
                cat("###############################################################\n")
                cat("Step ", i, ":\n")
                cat("stability values rho:\n")
                print(sort(rho[!is.na(rho)]))
                cat("variable with highest stability (smallest rho value):\t", Symbols[b[i]], "\n")
            }
        }
        names(R) <- 1:minNrHKs
        names(rho.min) <- 1:minNrHKs
        return(list(ranking = R, rho = rho.min))
    }else{
        stop("specified method not yet implemented")
    }
}
