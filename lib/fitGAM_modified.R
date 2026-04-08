fitGAM_modified <- function (counts, sds = NULL, pseudotime = NULL, 
        cellWeights = NULL, conditions = NULL, U = NULL, genes = seq_len(nrow(counts)), 
        weights = NULL, offset = NULL, nknots = 6, verbose = TRUE, 
        parallel = FALSE, BPPARAM = BiocParallel::bpparam(), 
        control = mgcv::gam.control(), sce = TRUE, family = "nb", 
        gcv = FALSE) 
    {
        if (is.null(counts)) 
            stop("Provide expression counts using counts", " argument.")
        if (is.null(sds) & (is.null(pseudotime) | is.null(cellWeights))) {
            stop("Either provide the slingshot object using the sds ", 
                "argument, or provide pseudotime and cell-level weights ", 
                "manually using pseudotime and cellWeights arguments.")
        }
        if (!is.null(sds)) {
            if (is(sds, "SlingshotDataSet") | is(sds, "PseudotimeOrdering")) {
                if (!sce) {
                  warning(paste0("If an sds argument is provided, the sce argument is ", 
                    "forced to TRUE "))
                  sce <- TRUE
                }
            }
            else stop("sds argument must be a SlingshotDataSet or ", 
                "PseudotimeOrdering object.")
            pseudotime <- slingPseudotime(sds, na = FALSE)
            cellWeights <- slingCurveWeights(sds)
        }
        if (any(is.na(pseudotime))) {
            stop("The pseudotimes contain NA values, and these cannot be used", 
                " for GAM fitting.")
        }
        if (any(is.na(cellWeights))) {
            stop("The cellWeights contain NA values, and these cannot be used", 
                " for GAM fitting.")
        }
        if (!is.null(conditions)) {
            if (class(conditions) != "factor") 
                stop("conditions must be a factor vector.")
            if (length(conditions) != ncol(counts)) {
                stop("conditions vector must have same length as number of cells.")
            }
            if (nlevels(conditions) == 1) {
                message("Only one condition was provided. Will run fitGAM without conditions")
                conditions <- NULL
            }
            if (!sce) {
                warning(paste0("If conditions, tradeSeq will return", 
                  " a SingleCellExperiment object"))
            }
        }
        gamOutput <- .fitGAM_modified(counts = counts, U = U, pseudotime = pseudotime, 
            cellWeights = cellWeights, conditions = conditions, 
            genes = genes, weights = weights, offset = offset, 
            nknots = nknots, verbose = verbose, parallel = parallel, 
            BPPARAM = BPPARAM, control = control, sce = sce, 
            family = family, gcv = gcv)
        if (!sce) {
            return(gamOutput)
        }
        sc <- SingleCellExperiment(assays = list(counts = counts[genes, 
            ]))
        SummarizedExperiment::colData(sc)$crv <- S4Vectors::DataFrame(pseudotime = pseudotime, 
            cellWeights = cellWeights)
        df <- tibble::enframe(gamOutput$Sigma, value = "Sigma")
        df$beta <- tibble::tibble(beta = gamOutput$beta)
        df$converged <- gamOutput$converged
        suppressWarnings(rownames(df) <- rownames(counts[genes, 
            ]))
        SummarizedExperiment::rowData(sc)$tradeSeq <- df
        if (is.null(conditions)) {
            SummarizedExperiment::colData(sc)$tradeSeq <- tibble::tibble(X = gamOutput$X, 
                dm = gamOutput$dm)
        }
        else {
            SummarizedExperiment::colData(sc)$tradeSeq <- tibble::tibble(X = gamOutput$X, 
                dm = gamOutput$dm, conditions = conditions)
        }
        S4Vectors::metadata(sc)$tradeSeq <- list(knots = gamOutput$knotPoints)
        return(sc)
    }

.checks <- function(pseudotime, cellWeights, U, counts, conditions, family) 
{
    if(family == "nb"){
        if (any(counts < 0)) {
      stop("All values of the count matrix should be non-negative")
    }
  }
    if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
        if (!identical(dim(pseudotime), dim(cellWeights))) {
            stop("pseudotime and cellWeights must have identical dimensions.")
        }
    }
    if (!is.null(U)) {
        if (!(nrow(U) == ncol(counts))) {
            stop("The dimensions of U do not match those of counts.")
        }
    }
    if (!is.null(dim(pseudotime)) & !is.null(dim(cellWeights))) {
        if (!identical(nrow(pseudotime), ncol(counts))) {
            stop("pseudotime and count matrix must have equal number of cells.")
        }
        if (!identical(nrow(cellWeights), ncol(counts))) {
            stop("cellWeights and count matrix must have equal number of cells.")
        }
    }
    if (!is.null(conditions)) {
        if (!is(conditions, "factor")) {
            stop("conditions must be a vector of class factor.")
        }
    }
    if (any(is.na(pseudotime)[cellWeights > 0])) {
        stop("Pseudotime contains NA values for non-zero weights.")
    }
    if (any(is.na(pseudotime))) {
        warning("Pseudotime contains NA values.")
    }
}

.assignCells <- function (cellWeights) 
{
    if (is.null(dim(cellWeights))) {
        if (any(cellWeights == 0)) {
            stop("Some cells have no positive cell weights.")
        }
        else {
            return(matrix(1, nrow = length(cellWeights), ncol = 1))
        }
    }
    else {
        if (any(rowSums(cellWeights) == 0)) {
            stop("Some cells have no positive cell weights.")
        }
        else {
            normWeights <- sweep(cellWeights, 1, FUN = "/", STATS = apply(cellWeights, 
                1, sum))
            wSamp <- apply(normWeights, 1, function(prob) {
                stats::rmultinom(n = 1, prob = prob, size = 1)
            })
            if (is.null(dim(wSamp))) {
                wSamp <- matrix(wSamp, ncol = 1)
            }
            else {
                wSamp <- t(wSamp)
            }
            return(wSamp)
        }
    }
}


.get_offset <- function (offset, counts) 
{
    if (is.null(offset)) {
        nf <- try(edgeR::calcNormFactors(counts), silent = TRUE)
        if (is(nf, "try-error")) {
            message("TMM normalization failed. Will use unnormalized library sizes", 
                "as offset.\n")
            nf <- rep(1, ncol(counts))
        }
        libSize <- colSums(as.matrix(counts)) * nf
        offset <- log(libSize)
        if (any(libSize == 0)) {
            message("Some library sizes are zero. Offsetting these to 1.\n")
            offset[libSize == 0] <- 0
        }
    }
    return(offset)
}

.findKnots <- function (nknots, pseudotime, wSamp) 
{
    for (ii in seq_len(ncol(pseudotime))) {
        assign(paste0("t", ii), pseudotime[, ii])
    }
    for (ii in seq_len(ncol(pseudotime))) {
        assign(paste0("l", ii), 1 * (wSamp[, ii] == 1))
    }
    tAll <- c()
    for (ii in seq_len(nrow(pseudotime))) {
        tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii, 
            ]))]
    }
    knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1))/(nknots - 
        1))
    if (any(duplicated(knotLocs))) {
        knotLocs <- stats::quantile(t1[l1 == 1], probs = (0:(nknots - 
            1))/(nknots - 1))
        if (any(duplicated(knotLocs))) {
            dupId <- duplicated(knotLocs)
            if (max(which(dupId)) == length(knotLocs)) {
                dupId <- duplicated(knotLocs, fromLast = TRUE)
                knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 
                  1], knotLocs[which(dupId) + 1]))
            }
            else {
                knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 
                  1], knotLocs[which(dupId) + 1]))
            }
        }
        if (any(duplicated(knotLocs))) {
            knotLocs <- seq(min(tAll), max(tAll), length = nknots)
        }
    }
    maxT <- max(pseudotime[, 1])
    if (ncol(pseudotime) > 1) {
        maxT <- c()
        for (jj in 2:ncol(pseudotime)) {
            maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l", 
                jj)) == 1])
        }
    }
    if (all(maxT %in% knotLocs)) {
        knots <- knotLocs
    }
    else {
        maxT <- maxT[!maxT %in% knotLocs]
        replaceId <- vapply(maxT, function(ll) {
            which.min(abs(ll - knotLocs))
        }, FUN.VALUE = 1)
        knotLocs[replaceId] <- maxT
        if (!all(maxT %in% knotLocs)) {
            warning(paste0("Impossible to place a knot at all endpoints.", 
                "Increase the number of knots to avoid this issue."))
        }
        knots <- knotLocs
    }
    knots[1] <- min(tAll)
    knots[nknots] <- max(tAll)
    knotList <- lapply(seq_len(ncol(pseudotime)), function(i) {
        knots
    })
    names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))
    return(knotList)
}