`nlt2` <-
function (x, f, J, Pred = 4, neighbours = 1, closest = FALSE, 
    intercept = TRUE, nkeep = 2, trule = "median",verbose=TRUE) 
{
    n <- length(x)
    vec <- matrix(0, J, n - 2)
    denf <- list()
    ghatmat <- matrix(0, J, n)
    ghatnat <- matrix(0, 1, n)
    for (i in 1:J) {
	if(verbose){
		cat(i,"...\n")
	}
        v <- sample(1:n, (n - 2), FALSE)
        vec[i, ] <- as.row(v)
        denf[[i]] <- denoiseperm2(x, f, pred = Pred, neigh = neighbours, 
            int = intercept, clo = closest, keep = nkeep, rule = trule, 
            per = v)
        ghatmat[i, ] <- as.row(denf[[i]]$fhat$coeff)
    }
    aveghat <- apply(ghatmat, 2, mean)
    df <- denoise2(x, f, pred = Pred, neigh = neighbours, int = 
intercept, 
        clo = closest, keep = nkeep, rule = trule)
    ghatnat <- as.row(df$fhat$coeff)
    return(list(vec = vec, ghatmat = ghatmat, ghatnat = ghatnat, 
        aveghat = aveghat))
}

