\name{Cens.SMN.PCR}
\alias{Cens.SMN.PCR}

\title{
Estimation and diagnostic for partially censored regression models
}
\description{
Return the MPL estimates obtained through a EM-type algorithm  for partially censored regression models under scale-mixture of normal (SMN) distributions (see, some members are the normal, student's-t, slash and contaminated normal distribution). The types of censoring considered are left and right. Graphics for diagnostic analysis such as case-deletion and local influence techniques are provided to show its robust aspect against outlying and influential observations.

}
\usage{
Cens.SMN.PCR(x, y, c, cens = "left", tt, nu = NULL, error = 10^-6, iter.max = 200,
type = "Normal", lambda.FIX = TRUE, nu.FIX = TRUE, lambda.in = 10^-3, k = 1,
Diagnostic = TRUE, a = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix or vector of covariates.
}
  \item{y}{
Vector of responses.
}
  \item{c}{
Vector of censoring indicators. For each observation: 1 if
censored and 0 if non-censored.
}
  \item{cens}{
'left' for left censoring and 'right' for rigth censoring.
}
  \item{tt}{
Vector of values of a continuous covariate for the non-parametric component of the model.
}
  \item{nu}{
Initial value of the parameter of the SMN family. In the case of the student's-t and slash is a scalar, in the contaminated normal is a vector bidimensional.
}
  \item{error}{
The convergence maximum error. By default = 10^-6.
}
  \item{iter.max}{
The maximum number of iterations of the EM-type algorithm. By default = 200.
}
  \item{type}{
Represents the type of distribution to be used in fitting: 'Normal' for normal, 'T' for student's-t, 'Slash' for slash and 'NormalC' for contaminated normal distribution respectively. By default ='Normal'
}
  \item{lambda.FIX}{
\code{TRUE} or \code{FALSE}. Indicate if smoothing parameter  will be estimated. By default = \code{TRUE}.
}
  \item{nu.FIX}{
\code{TRUE} or \code{FALSE}. Indicate if \eqn{\nu}  will be estimated. By default = \code{TRUE}.
}
  \item{lambda.in}{
Initial value of smoothing parameter (also known in the literature as \eqn{\alpha}).
}
  \item{k}{
For the local influence in explanatory variable perturbation, indicates the \eqn{k-}th explanatory variable (assumed continuous) of the design matrix  \eqn{X} to be perturbed.
}
  \item{Diagnostic}{
\code{TRUE} or \code{FALSE}. Indicates if diagnostic graph should be built for the fitted model (index plot in local influence). By default = \code{TRUE}.
}
  \item{a}{
The value for \eqn{a} considered in the benchmark value for the index plot in local influence: \eqn{M(0)_{l} > \bar{M(0)}+a*SM(0)}.
}
}
\details{
We consider a partial linear model which belongs to the class of semiparametric regression with vector of response  \eqn{Y=(Y_{1},...,Y_{n})} and with errors \eqn{\epsilon_{i}} which are independent and identically distributed according to a SMN distribution. To be more precise,

\deqn{Y_{i} = x_i^{T}\beta +n_i^{T}f + \epsilon_{i},}

for \eqn{i=1,...,n}, where \eqn{f = (f(t_1^{0}),...,f(t_r^{n})^{T}} is an \eqn{r x 1} vector with \eqn{t_1^{0},...,t_r^{n}} being the distinct and ordered values of \eqn{t_i}; \eqn{n_i} is a \eqn{r x 1} vector of incidence whose \eqn{s}-th element equals the indicator function \eqn{I(t_i=t_s^{0})} for \eqn{s=1,...,r}.
}
\value{

\item{beta}{EM estimates for the parametric component.}
\item{sigma2}{EM estimates for the scale parameter.}
\item{lambda}{If \code{lambda.FIX} = \code{FALSE}, it returns the estimated value of the smoothing parameter, else returns the initial value assigned in \code{lambda.in}.}
\item{AIC}{AIC criteria for model selection.}
\item{ff}{EM estimates for the nonparametric component.}
\item{yest}{Predicted values of the model.}
\item{loglik}{Value of the log-likelihood under the fitted model.}
\item{iter}{Number of iterations of the EM algorithm.}
\item{nu}{If \code{nu.FIX} = \code{FALSE}, it returns the estimated value of \eqn{\nu} parameter, else returns the initial value assigned in \code{nu}.}
\item{MI}{Observed information matrix.}
\item{D}{A list of objects for diagnostic analysis that contains: the Hessian matrix (\code{Hessian}), values for generalized Cook's distance (\code{GD1}) and the values of the conformal normal curvature for the following perturbation schemes: Case-weight (\code{Curvature_W}), scale (\code{Curvature_S}), explanatory variable (\code{Curvature_E}) and response variable (\code{Curvature_R}).}

}
\references{
Ferreira, C. S., & Paula, G. A. (2017). Estimation and diagnostic for skew-normal partially linear models. Journal of Applied Statistics, 44(16), 3033-3053.

Ibacache-Pulgar, G., Paula, G. A., & Cysneiros, F. J. A. (2013). Semiparametric additive models under symmetric distributions. Test, 22(1), 103-121.

Ibacache-Pulgar, G., & Paula, G. A. (2011). Local influence for Student-t partially linear models. Computational Statistics & Data Analysis, 55(3), 1462-1478.
}
\author{
Marcela Nunez Lemus, Christian E. Galarza, Larissa Avila Matos and Victor H. Lachos
}

\section{Warning}{For the contaminated normal case, if nu parameters were close to the bounds, i.e., close to 0 or 1, computational problems could arrise.}

\note{When \code{lambda.FIX = FALSE} the algorithm may take a long time to converge.}

\seealso{
\code{\link[SMNCensReg]{CensReg.SMN}}
}
\examples{
dtawage = get(data(PSID1976,package = "AER"))
y  = dtawage$wage
cc = c(rep(0,428),rep(1,325))
tt = dtawage$exper
x  = cbind(dtawage$education,dtawage$age, dtawage$hhours, dtawage$hwage, dtawage$tax,
dtawage$youngkids, dtawage$oldkids)

#Normal case by default with only 10 iterations
PCR.default1 = Cens.SMN.PCR(x=x, y=y, c=cc, cens="left",tt =tt,iter.max = 10,Diagnostic = FALSE)

\dontrun{
#This may take few minutes
#Normal case by default with full (200) iterations
PCR.default2 = Cens.SMN.PCR(x=x, y=y, c=cc, cens="left",tt =tt)

#contaminated normal case
PCR.CN = Cens.SMN.PCR(x=x, y=y, c=cc, cens="left",tt =tt,type="NormalC",
nu = c(0.1,0.1),iter.max = 100)
}
}
\keyword{package}
\keyword{censored}
\keyword{censoring}
\keyword{diagnostic}
\keyword{partial}
\keyword{partially}
