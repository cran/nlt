transmatdualperm2 <-
function (x, f, Pred = 5, neigh = 1, int = TRUE, clo = TRUE, 
    keep = 2, perm=sample(1:length(x),(length(x)-keep),FALSE), varonly = 
    FALSE) 
{
    Wnew <- NULL
    if (length(x) == length(unique(x))) {
        out <- fwtnpperm2(x, f, LocalPred = Pred, neighbours = neigh, 
            intercept = TRUE, closest = clo, nkeep = keep,mod=perm)
    }
    else {
        out <- fwtnpmp(x, f, LocalPred = Pred, neighbours = neigh, 
            intercept = TRUE, closest = clo, nkeep = keep, mpdet = "min")
    }
    n <- length(out$X)
    Alistdual <- list()
    Tlistdual <- list()
    matno <- n - keep
    lca <- out$lca
    pointsin <- out$pointsin
    remlist <- lca[, 1]
    newpoints <- (c(pointsin, rev(remlist)))
    if (matno > 0) {
    		nn<-lca[matno,2]
	    Alistdual[[1]] <- Amatdual2(matno, pointsin, remlist,lca[matno,3:(2+nn)], lca[matno,(3+2*nn):(2+3*nn)],lca[matno,(3+nn):(2+2*nn)])
            Tlistdual[[1]] <- Alistdual[[1]]
            if (matno > 1) {
	            for (i in 2:matno) {
            		nn <- lca[matno-i+1, 2]
            		nbrs <- lca[matno-i+1, 3:(2 + nn)]
            		alpha <- lca[matno-i+1, (3 + nn):(2 + 2 * nn)]
            		weights <- lca[matno-i+1, (3 + 2 * nn):(2 + 3 * nn)]
            		Ai <- Amatdual2(matno-i+1, pointsin, remlist, nbrs, weights, alpha)
            		Alistdual[[i]] <- Ai
      		        augment <- rbind(cbind(Tlistdual[[i - 1]], 0),  0)
      			augment[nrow(augment), nrow(augment)] <- 1
      			Tlistdual[[i]] <- augment %*% Alistdual[[i]]
        	    }
      	    }

        W <- Tlistdual[[matno]]
        x <- out$X
        Wnew <- matrix(0, length(x), length(x))
        reo <- NULL
        reo <- match(x, x[newpoints])
         Wnew<-W[reo,reo]
    }
    if (varonly) {
        return(diag(Wnew))
    }
    else {
        return(list(out = out, Wnew = Wnew, x = x))
    }
}

