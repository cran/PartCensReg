Cens.SMN.PCR <-
function(x, y, c, cens="left", tt, nu=NULL, error=10^-6, iter.max=200, type = "Normal", lambda.FIX = TRUE, nu.FIX = TRUE, lambda.in = 10^-3, k = 1, Diagnostic = TRUE, a = 2)
{
  #== Required packages

  # require("ssym")
  # require("optimx")
  # require("MASS")
  # require("Matrix")

  #---------------------------------------------------------------------#
  #                              Validations                            #
  #---------------------------------------------------------------------#


  #== Matrix column labels

  namesx = ('x1     ')
  if(ncol(as.matrix(x))>1)
  {
    for(i in 2:ncol(as.matrix(x))){namesx = cbind(namesx, paste("x",i,"     ",sep=""))}
  }

  #== Validation, dimension of the dataset

  if(ncol(as.matrix(y)) > 1) stop("Only univariate partially regression models, argument must be one-dimensional.")
  if(ncol(as.matrix(tt)) > 1) stop("Only univariate partially regression models, argument must be one-dimensional.")
  if(ncol(as.matrix(c)) > 1) stop("Only univariate partially regression models, argument must be one-dimensional.")
  if( length(y) != nrow(as.matrix(x)) ) stop("The number of rows in the matrix X it must be the same than Y.")
  if( length(c) != nrow(as.matrix(x)) ) stop("The number of rows in the matrix X it must be the same than c.")
  if( length(tt) != nrow(as.matrix(x)) ) stop("The number of rows in the matrix X it must be the same than tt.")

  #== Validation, no data and NA's

  if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
  if(sum(is.na(y)==TRUE) > 0) stop("NA's values in y")
  if(sum(is.na(tt)==TRUE) > 0) stop("NA's values in tt")
  if(sum(is.na(x)==TRUE) > 0) stop("NA's values in X")

  #== Validation, type of distribution

  if( (type != "Normal") && (type != "T") && (type != "Slash") && (type != "NormalC")) stop("Family not recognized. Please check documentation.")

  #== Validation, type of censoring

  if( (cens != "left") && (cens != "right"))  stop("Censoring not recognized: 'left' for left censoring and 'right' for right censoring.")

  if(cens=="left"){cens = 1}else{cens = 2}

  #== Validation, parameters of the distribution belonging to the SMN family

  if(type == "Normal" && !is.null(nu)){warning("Nu parameter not considered for normal case.",immediate. = TRUE)}


  if( type == "T")
  {
    if(length(nu) > 1 | nu <= 0) stop("For the Student's-t distribution, nu parameter must be a positive scalar")
    if(is.null(nu)) stop("nu parameter must be provided.")
  }

  if( type == "Slash")
  {
    if(length(nu) > 1 | nu <= 0) stop("For the Slash distribution, nu parameter must be a positive scalar")
    if(is.null(nu)) stop("nu parameter must be provided.")
  }

  if(type == "NormalC")
  {
    warning("When nu parameters are close to the bounds, i.e., 0 or 1, computational problems could arrise.",immediate. = TRUE)
    if(length(nu) != 2) stop("For the Contaminated Normal distribution, nu must be a bidimensional vector.")
    if(nu[1] <=0 | nu[1] >= 1) stop("nu1 must belong to the interval (0,1)")
    if(nu[2] <=0 | nu[2] >= 1) stop("nu2 must belong to the interval (0,1)")
  }

  #== Validation, smoothing parameter lambda > 0

  if(lambda.in <= 0) stop("lambda parameter must be positive.")

  #== Validation, additional parameters

  if(!any(k == seq(1,ncol(as.matrix(x))))) stop("k must be positive integer <= ncol(x).")
  if(a <= 0) stop("a must be positive constant.")

  #== Validation, arguments support of the function

  if(iter.max <= 0 | iter.max%%1 != 0) stop("iter.max must be a positive integer.")
  if(error <=0 | error > 1) stop("error must belong to the interval (0,1]")

  if(!is.logical(lambda.FIX) | !is.logical(nu.FIX) | !is.logical(Diagnostic)) stop("Parameters lambda.FIX, nu.FIX and Diagnostic must be logical (TRUE/FALSE) variables.")

  #---------------------------------------------------------------------#
  #                              EM outputs                             #
  #---------------------------------------------------------------------#


  out.EM = EMSpline.censFinal_MobsSMN(x, y, c, cens, tt, nu, error, iter.max, type, delta.in=NA, lambda.FIX, nu.FIX, lambda.in, k)

  betas    = round(out.EM$beta, 4)
  sigma2   = round(out.EM$sigma2, 4)
  MI_obs   = sqrt(diag(solve(out.EM$MI)))
  SEbeta   = MI_obs[1:length(betas)]
  SEsigma2 = MI_obs[length(betas)+length(out.EM$ff)+1]
  SE       = round(c(SEbeta,SEsigma2), 4)

  Estimates      = cbind(rbind(betas,sigma2),SE)
  namesEstimates = colnames(x)
  colx           = ncol(as.matrix(x))

  greeks = c(alpha='\u03b1', sigma='\u03c3\u00B2', nu='\u03BD')

  if(length(namesEstimates)==0) namesEstimates = namesx[1:colx]
  dimnames(Estimates) = list(c(namesEstimates,paste0(greeks['sigma'])),c("Estimates", "SE"))


  if( (type=="T") || (type=="Slash"))
  {
    param1            = matrix(round(out.EM$nu, 5),ncol=1,nrow=1)
    dimnames(param1)  = list(c(paste0(greeks['nu'])),"")
  }

  if( type=="NormalC")
  {
    param2            = matrix(round(out.EM$nu, 5),ncol=1,nrow=2)
    dimnames(param2)  = list(c(paste0(greeks['nu'],"1"), paste0(greeks['nu'],"2")),"")
  }

 if( type=="Normal")
  {
    alpha            = t(as.matrix(out.EM$lambda))
    row.names(alpha) = paste0(greeks['alpha'])
    colnames(alpha)  = " "
  }

  if( (type=="T") || (type=="Slash"))
  {
    alpha1            = matrix(out.EM$lambda,ncol=1,nrow=1)
    dimnames(alpha1)  = list(c(paste0(greeks['alpha'])),"")
  }

 if( type=="NormalC")
  {
    alpha2            = matrix(out.EM$lambda,ncol=1,nrow=1)
    dimnames(alpha2)  = list(c(paste0(greeks['alpha'])),"")
  }

  cat('\n')
  cat('--------------------------------------------------------------\n')
  cat('     Partially censored regression models with SMN errors     \n')
  cat('--------------------------------------------------------------\n')
  print(Estimates)

 if(type!="Normal")
  {
    if(type=="T"|type=="Slash")
    {
      print(param1)
    }
    else
    {
      print(param2)
    }
  }

   if ( type=="Normal") {
   print(alpha)
     } else if ( type=="T"|type=="Slash") {
   print(alpha1)
     } else
   print(alpha2)

   cat('--------------------------------------------------------------\n')
   cat('\r \n')
   criteriaPCR = c(out.EM$loglik, out.EM$AIC)
   criteriaFin = round(t(as.matrix(criteriaPCR)),digits=3)
   dimnames(criteriaFin) = list(c("Value"),c("Loglik", "AIC"))
   cat('\n')
   cat('Model selection criteria\n')
   cat('------------------------------------\n')
   print(criteriaFin)
   cat('------------------------------------\n')
   cat('\r \n')


  if(Diagnostic=="TRUE")
    {
     CaseDeletion      = plot(out.EM$D$GD1, las=1, ylab=expression(paste("GD"[i]^{1})))
     dev.new()
     par(mfrow=c(2,2), mar=c(4,4.5,1.5,1.5) + 0.1)
     WeightScheme      = plot(out.EM$D$Curvature_W, las=1, ylab="Case weight, M(0)")
     abline(h=mean(out.EM$D$Curvature_W)+a*sd(out.EM$D$Curvature_W), col="red", lty=2)
     ScaleScheme       = plot(out.EM$D$Curvature_S, las=1, ylab="Case scale, M(0)")
     abline(h=mean(out.EM$D$Curvature_S)+a*sd(out.EM$D$Curvature_S), col="red", lty=2)
     ExplanatoryScheme = plot(out.EM$D$Curvature_E, las=1, ylab="Case explanatory, M(0)")
     abline(h=mean(out.EM$D$Curvature_E)+a*sd(out.EM$D$Curvature_E), col="red", lty=2)
     ResponseScheme    = plot(out.EM$D$Curvature_R, las=1, ylab="Case response, M(0)")
     abline(h=mean(out.EM$D$Curvature_R)+a*sd(out.EM$D$Curvature_R), col="red", lty=2)
    }

out.EM
}
