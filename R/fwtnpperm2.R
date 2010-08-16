fwtnpperm2 <-
function (input, f, nkeep = 2, intercept = TRUE, initboundhandl = 0, 
    neighbours = 1, closest = FALSE, LocalPred = 1, mod = 
    sample(1:length(input), (length(input) - 
            nkeep), FALSE),do.W = FALSE, 
    varonly = FALSE) 
{
    n <- length(f)
    if (do.W) {
        if (varonly) {
            v <- rep(0, times = n)
        }
        else {
            v <- 0
        }
        W <- matrix(0, n, n)
    }
    else {
        v <- 0
        W <- 0
    }
    if ((n - nkeep) > 0) {
        coeff <- rep(0, times = n)
        pointsin <- coeff
        lengthsremove <- rep(0, times = n - nkeep)
        lca <- matrix(0, n - nkeep, 6 * neighbours + 5)
        ans <- .C("fwtnp", as.double(input), as.double(f), as.integer(nkeep), 
            as.integer(intercept), as.integer(initboundhandl), 
            as.integer(neighbours), as.integer(closest), as.integer(LocalPred), 
            as.integer(n), coeff = as.double(coeff), lr = as.double(lengthsremove), 
            lengths = as.double(coeff), lca = as.double(t(lca)), 
            po = as.integer(pointsin), nc = as.integer(0), 
            mod=as.integer(mod),as.integer(do.W), 
            W = as.double(t(W)), as.integer(varonly), v = as.double(v), 
            PACKAGE = "adlift")
        ans$po <- ans$po[1:nkeep]
        ans$lengths <- ans$lengths[1:nkeep]
        ans$lca <- matrix(ans$lca[1:((n - nkeep) * ans$nc)], 
            ncol = ans$nc, byrow = TRUE)
        coeff <- ans$coeff
        lengthsremove <- ans$lr
        lca <- ans$lca
        pointsin <- ans$po
        lengths <- ans$lengths
        if(varonly){
        	Wnew<-matrix(ans$W[1:(nkeep*n)],byrow=TRUE,nrow=nkeep)
        }
        else{
	        Wnew <- matrix(ans$W, n, n, byrow = TRUE)
	}        
        v <- ans$v
    }
    else {
        coeff <- f
        pointsin <- order(input)
        lengths <- lengthsremove <- lca <- NULL
    }
    if (varonly) {
        Wnew <- NULL
    }
    return(list(X = input, coeff = coeff, lengthsremove = lengthsremove, 
        lca = lca, pointsin = pointsin, lengths = lengths, W = Wnew, 
        v = v))
}

