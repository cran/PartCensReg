\name{PartCensReg-package}
\alias{PartCensReg-package}
\docType{package}
\title{
\packageTitle{PartCensReg}
}
\description{
\packageDescription{PartCensReg}
}

\references{
Ferreira, C. S., & Paula, G. A. (2017). Estimation and diagnostic for skew-normal partially linear models. Journal of Applied Statistics, 44(16), 3033-3053.

Ibacache-Pulgar, G., Paula, G. A., & Cysneiros, F. J. A. (2013). Semiparametric additive models under symmetric distributions. Test, 22(1), 103-121.

Ibacache-Pulgar, G., & Paula, G. A. (2011). Local influence for Student-t partially linear models. Computational Statistics & Data Analysis, 55(3), 1462-1478.
}

\keyword{package}
\keyword{censored}
\keyword{censoring}
\keyword{diagnostic}
\keyword{partial}
\keyword{partially}

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
