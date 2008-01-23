\name{fwtnpperm}
\alias{fwtnpperm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform a modified version of the lifting algorithm}
\description{
  Given a signal, a chosen prediction strategy and neighbourhood size and configuration, this function applies the lifting algorithm that removes one coefficient at a time following a prespecified order of point removal (or path, trajectory), rather than letting the algorithm decide the removal order. We refer to this type of lifting algorithm as the modified lifting scheme.
}
\usage{
fwtnpperm(input, f, nkeep = 2, intercept = TRUE, initboundhandl = "reflect", updateboundhandl = "add", neighbours = 1, closest = FALSE, LocalPred = AdaptPred, mod = sample(1:length(input),(length(input)-nkeep),FALSE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{Vector of any length (possibly irregularly spaced) that gives the grid locations at which the signal is observed.}
  \item{f}{Vector of the same length as \code{input} that gives the signal values corresponding to the grid locations.}
  \item{nkeep}{Number of scaling points we want after processing the signal by means of the modified lifting scheme. The usual choice is \code{nkeep}=2.}
  \item{intercept}{Specifies whether (\code{intercept}=TRUE) or not (\code{intercept}=FALSE) an intercept is to be used in the regression curve. For \code{LocalPred}=\code{AdaptPred} or \code{AdaptNeigh}, the algorithm automatically makes this choice.}
  \item{initboundhandl}{Character string, it indicates the strategy used for handling the interval boundaries at the start of the transform: \code{"reflect"} indicates that the intervals corresponding to the first and last data points are taken to have the respective grid values as midpoints, while \code{"stop"} indicates that the first and last intervals have as outer endpoints the grid points themselves.}
  \item{updateboundhandl}{Character string, it indicates the strategy used for updating the intervals associated to the neighbours of the points removed from the boundary:\code{"reflect"} indicates that the neighbouring intervals are symmetrical around the corresponding grid point, \code{"stop"} indicates that the neighbouring intervals end at the respective data points, and \code{"add"} indicates that the neighbouring intervals have outer endpoints where the boundary intervals did.}
  \item{neighbours}{Number of neighbours to be used for defining the neighbourhood of each point that has to be removed. If (\code{closest}=FALSE), then this gives the number of neighbours to be used on each side of the removed point.}
  \item{closest}{If (\code{closest}=TRUE) or (\code{closest}=FALSE), then at each step the neighbours are in closest, respectively symmetrical configuration.}
  \item{LocalPred}{The type of regression to be used in the prediction step of the modified lifting algorithm. Choices are linear, quadratic or cubic (respectively, \code{LinearPred}, \code{QuadPred} or \code{CubicPred}), or two adaptive procedures which automatically choose the degree used in regression, (\code{AdaptPred} or \code{AdaptNeigh}).}
  \item{mod}{Vector of length (length(\code{x})-\code{nkeep}) that gives the order of point removal to be used by the modified lifting algorithm. This is in fact the trajectory or path previously mentioned.}
}
\details{
 The function here is not designed to cope with multiple observations at the same grid point, although it can be modified to do so. If multiple data is supplied, then the algorithm performs the usual lifting algorithm coded in the function \code{fwtnpmp} (designed to handle multiple entries) of the package 'adlift'. The algorithm can be easily inverted by using the function \code{invtnp} of the package 'adlift'. Also the functions that appear as options for \code{LocalPred} can be found in the same 'adlift' package.
}
\value{
  \item{X}{Vector of grid points used in the transform.}
  \item{coeff }{Vector (in the order of \code{X}) giving the detail and scaling coefficients obtained after lifting the initial signal.}
  \item{origlengths}{Vector of initial interval lengths corresponding to each gridpoint.}
  \item{lengths}{Vector of updated interval lengths corresponding to the scaling coefficients. Its length is \code{nkeep}.}
  \item{lengthsremove}{Vector giving the (updated) interval lengths corresponding to each removed point at its removal time. It corresponds to the removal order given by \code{removelist}.}
  \item{pointsin}{Vector of indices into \code{X} giving the scaling coefficients after the lifting algorithm removed the points specified by the trajectory. The length of this vector is \code{nkeep}.}
  \item{removelist}{Vector of indices into \code{X} giving the lifted coefficients in their order of removal. This should coincide with the given trajectory.}
  \item{neighbrs}{List of indices into \code{X}. Each list entry gives the indices of the neighbours of the removed point used at that particular step of the transform.}
  \item{neighbours}{Number of neighbours used in the prediction step of the transform.}
  \item{gamlist}{List of prediction weights used in the prediction step of the transform.}
  \item{alphalist}{List of update coefficients used in the update step of the decomposition.}
  \item{schemehist}{Vector of character strings indicating the type of regression used at each step of the transform. This is NULL when a non-adaptive regression step is chosen.}
  \item{interhist}{Boolean vector indicating whether or not an intercept was used in the regression curve at each step. This is NULL when a non-adaptive regression step is chosen.}
  \item{clolist}{Boolean vector recording the neighbour configuration chosen at each step: TRUE corresponds to closest neighbours, while FALSE to symmetrical neighbours. It is NULL except for \code{LocalPred}=\code{AdaptNeigh}.}
}
\references{See the paper 'A "nondecimated" lifting transform.' by Knight, M.I. and Nason, G.P. (2008) for further details.}
\author{Marina Knight (marina.knight@bristol.ac.uk)}
\note{Use this function together with the "adlift" package available from CRAN.}
\seealso{\code{fwtnp},\code{fwtnpmp},\code{invtnp} of package 'adlift'}
\examples{
#construct an irregular grid of length 256
input<-runif(256)

#now construct a Doppler signal on the grid above 
f<-make.signal2("doppler",x=input)

#generate the vector which will indicate the order of point removal that will be followed by the modified lifting algorithm (i.e. the trajectory)
#vec below gives the first (length(input)-nkeep) entries of a random permutation of (1:length(input)) 
vec<-sample(1:256,254,FALSE)

#now transform the initial signal (input,f) into a set of wavelet and scaling coefficients by using a modified lifting transform 
#that follows the removal order in vec and linear regression prediction steps and neighbourhoods of size 2 in symmetrical configuration
#the set of scaling and detail coefficients can be found in out$coeff
out<-fwtnpperm(input,f,LocalPred=LinearPred,neighbours=1,closest=FALSE,intercept=TRUE,nkeep=2,mod=vec) 
}
\keyword{methods}
