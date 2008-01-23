\name{transmatdualperm}
\alias{transmatdualperm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the matrix associated to the modified lifting algorithm}
\description{
 Given a signal, a method of prediction (including neighbourhood definition) and a trajectory for the modified lifting algorithm to follow, this function computes the matrix that represents the transform from the raw data to the set of detail and scaling coefficients.   
}
\usage{
transmatdualperm(x, f, Pred = AdaptPred, neigh = 1, int = TRUE, clo = FALSE, keep = 2, perm = sample(1:length(x),(length(x)-keep),FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of any length (possibly irregularly spaced) that gives the grid locations at which the signal is observed.}
  \item{f}{Vector of the same length as \code{x} that gives the signal values corresponding to the \code{x}-locations.}
  \item{Pred}{The type of regression to be used in the prediction step of the modified lifting algorithm. Choices are linear, quadratic or cubic (respectively, \code{LinearPred}, \code{QuadPred} or \code{CubicPred}), or two adaptive procedures which automatically choose the degree used in regression, (\code{AdaptPred} or \code{AdaptNeigh}).}
  \item{neigh}{Number of neighbours to be used for defining the neighbourhood of each point that has to be removed. If (\code{clo}=FALSE), then this gives the number of neighbours to be used on each side of the removed point.}
  \item{int}{Specifies whether (\code{int}=TRUE) or not (\code{int}=FALSE) an intercept is to be used in the regression curve. For \code{Pred}=\code{AdaptPred} or \code{AdaptNeigh}, the algorithm automatically makes this choice.}
  \item{clo}{If (\code{clo}=TRUE) or (\code{clo}=FALSE), then at each step the neighbours are in closest, respectively symmetrical configuration.}
  \item{keep}{Number of scaling points we want at the end of the application of the modified lifting scheme. The usual choice is \code{keep}=2.}
  \item{perm}{Vector of length (length(\code{x})-\code{keep}). It gives the trajectory for the modified lifting algorithm to follow, i.e. it gives the order of point removal.}
}
\details{
  The matrix associated to the modified lifting transform is of size (length(\code{x}),length(\code{x})), and its columns correspond to the points in their order observation. The functions that appear as choices of \code{Pred} can be found in the package 'adlift'.
}
\value{
  \item{out}{Outcome of the modified lifting algorithm (\code{\link{fwtnpperm}}).}
  \item{Wnew}{Matrix associated to the modified lifting transform.}
  \item{x}{Vector of grid locations, the same as the \code{x} in \code{Arguments}.}
}
\references{See the paper 'A "nondecimated" lifting transform.' by Knight, M.I. and Nason, G.P. (2008) for further details.}
\author{Marina Knight (marina.knight@bristol.ac.uk)}
\note{Use this function together with the "adlift" package available from CRAN.}
\seealso{\code{\link{fwtnpperm}}}
\examples{
#construct the grid locations
x<-runif(256)

#construct a signal 
f<-runif(256) 

#generate the vector which will indicate the order of point removal that will be followed by the modified lifting algorithm (i.e. the trajectory)
#vec below gives the first (length(x)-keep) entries of a random permutation of (1:length(x)) 
vec<-sample(1:256,254,FALSE)  

#now you can transform the initial signal (x,f) into a set of wavelet and scaling coefficients, and work out the matrix associated to the transform 
#in the example below, the lifting transform follows the path given by vec and adaptive prediction steps with neighbourhoods of sizes at most 2 
#the matrix associated to the transform is given by out$Wnew, and the output from the modified lifting transform is found in out$out
out<-transmatdualperm(x,f,Pred=AdaptNeigh,neigh=1,clo=TRUE,keep=2,perm=vec) 
}
\keyword{methods}
