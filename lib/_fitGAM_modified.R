.fitGAM_modified <- function (counts, U = NULL, pseudotime, cellWeights, conditions, 
    genes = seq_len(nrow(counts)), weights = NULL, offset = NULL, 
    nknots = 6, verbose = TRUE, parallel = FALSE, BPPARAM = BiocParallel::bpparam(), 
    aic = FALSE, control = mgcv::gam.control(), sce = TRUE, family = "nb", 
    gcv = FALSE) 
{
    if (!is.null(conditions)) {
        message("Fitting lineages with multiple conditions. This method has ", 
            "been tested on a couple of datasets, but is still in an ", 
            "experimental phase.")
    }
    if (is(genes, "character")) {
        if (!all(genes %in% rownames(counts))) {
            stop("The genes ID is not present in the models object.")
        }
        if (any(duplicated(genes))) {
            stop("The genes vector contains duplicates.")
        }
        id <- match(genes, rownames(counts))
    }
    else {
        id <- genes
    }
    if (parallel) {
        BiocParallel::register(BPPARAM)
        if (verbose) {
            BPPARAM$tasks = as.integer(40)
            BPPARAM$progressbar = TRUE
        }
    }
    if (is.null(dim(pseudotime))) {
        pseudotime <- matrix(pseudotime, nrow = length(pseudotime))
    }
    if (is.null(dim(cellWeights))) {
        cellWeights <- matrix(cellWeights, nrow = length(cellWeights))
    }
    .checks(pseudotime, cellWeights, U, counts, conditions, family)
    wSamp <- .assignCells(cellWeights)
    for (ii in seq_len(ncol(pseudotime))) {
        assign(paste0("t", ii), pseudotime[, ii])
    }
    for (ii in seq_len(ncol(pseudotime))) {
        assign(paste0("l", ii), 1 * (wSamp[, ii] == 1))
    }
    offset <- .get_offset(offset, counts)
    if (is.null(U)) {
        U <- matrix(rep(1, nrow(pseudotime)), ncol = 1)
    }
    knotList <- .findKnots(nknots, pseudotime, wSamp)
    teller <- 0
    converged <- rep(TRUE, length(genes))
    counts_to_Gam <- function(y) {
        teller <<- teller + 1
        nknots <- nknots
        if (!is.null(weights)) 
            weights <- weights[teller, ]
        if (!is.null(dim(offset))) 
            offset <- offset[teller, ]
        if (is.null(conditions)) {
            smoothForm <- stats::as.formula(paste0("y ~ -1 + U + ", 
                paste(vapply(seq_len(ncol(pseudotime)), function(ii) {
                  paste0("s(t", ii, ", by=l", ii, ", bs='cr', id=1, k=nknots)")
                }, FUN.VALUE = "formula"), collapse = "+"), " + offset(offset)"))
        }
        else {
            for (jj in seq_len(ncol(pseudotime))) {
                for (kk in seq_len(nlevels(conditions))) {
                  lCurrent <- get(paste0("l", jj))
                  id1 <- which(lCurrent == 1)
                  lCurrent[id1] <- ifelse(conditions[id1] == 
                    levels(conditions)[kk], 1, 0)
                  assign(paste0("l", jj, "_", kk), lCurrent)
                }
            }
            smoothForm <- stats::as.formula(paste0("y ~ -1 + U + ", 
                paste(vapply(seq_len(ncol(pseudotime)), function(ii) {
                  paste(vapply(seq_len(nlevels(conditions)), 
                    function(kk) {
                      paste0("s(t", ii, ", by=l", ii, "_", kk, 
                        ", bs='cr', id=1, k=nknots)")
                    }, FUN.VALUE = "formula"), collapse = "+")
                }, FUN.VALUE = "formula"), collapse = "+"), " + offset(offset)"))
        }
        s <- mgcv::s
        m <- suppressWarnings(try(withCallingHandlers({
            mgcv::gam(smoothForm, family = family, knots = knotList, 
                weights = weights, control = control)
        }, error = function(e) {
            converged[teller] <<- FALSE
            return(structure("Fitting errored", class = c("try-error", 
                "character")))
        }, warning = function(w) {
            converged[teller] <<- FALSE
        }), silent = TRUE))
        return(m)
    }
    if (parallel) {
        #expr_list <- split(as.matrix(counts)[id, ], row(as.matrix(counts)[id, ]))  # list of gene expression vectors
        #gamList <- BiocParallel::bplapply(expr_list, counts_to_Gam, BPPARAM = BPPARAM)
        gamList <- BiocParallel::bplapply(as.list(as.data.frame(t(as.matrix(counts)[id, ]))), counts_to_Gam, BPPARAM = BPPARAM) #原本输入的是矩阵，但BatchtoolsParam需要list
    }

    
    else {
        if (verbose) {
            gamList <- pbapply::pblapply(as.data.frame(t(as.matrix(counts)[id, 
                ])), counts_to_Gam)
        }
        else {
            gamList <- lapply(as.data.frame(t(as.matrix(counts)[id, 
                ])), counts_to_Gam)
        }
    }
    if (aic) {
        aicVals <- unlist(lapply(gamList, function(x) {
            if (class(x)[1] == "try-error") 
                return(NA)
            x$aic
        }))
        if (gcv) {
            gcvVals <- unlist(lapply(gamList, function(x) {
                if (class(x)[1] == "try-error") 
                  return(NA)
                x$gcv.ubre
            }))
            return(list(aicVals, gcvVals))
        }
        else return(aicVals)
    }
    if (sce) {
        betaAll <- lapply(gamList, function(m) {
            if (is(m, "try-error")) {
                beta <- NA
            }
            else {
                beta <- matrix(stats::coef(m), ncol = 1)
                rownames(beta) <- names(stats::coef(m))
            }
            return(beta)
        })
        betaAllDf <- data.frame(t(do.call(cbind, betaAll)))
        rownames(betaAllDf) <- rownames(counts)[id]
        SigmaAll <- lapply(gamList, function(m) {
            if (is(m, "try-error")) {
                Sigma <- NA
            }
            else {
                Sigma <- m$Vp
            }
            return(Sigma)
        })
        element <- min(which(!is.na(SigmaAll)))
        m <- gamList[[element]]
        X <- stats::predict(m, type = "lpmatrix")
        dm <- m$model[, -1]
        knotPoints <- m$smooth[[1]]$xp
        return(list(beta = betaAllDf, Sigma = SigmaAll, X = X, 
            dm = dm, knotPoints = knotPoints, converged = converged))
    }
    else {
        return(gamList)
    }
}

