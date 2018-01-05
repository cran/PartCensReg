dignostic.EM_Sem = function(x, y, tt, beta, ff, sigma2, lambda,  u0, u1, u2, k , type){

  p  = ncol(x)
  n  = nrow(x)

  tr = sort(unique(tt))
  r  = length(tr)
  t2 = ncs(tr, nknots=r)
  K  = attr(t2, "K")
  N  = matrix(data = 0,nrow = n,ncol = r)
  for(i in 1:r){N[tt==tr[i],i] = 1}
  qq = ncol(N)

  alpha = as.vector(c(beta,ff))
  theta = c(beta,ff,sigma2)

  soma1.1 = matrix(0,p,p)                   # 2a derived of beta
  soma2.2 = matrix(0,qq,qq)                 # 2a derived of ff
  soma3.3 = 0                               # 2a derived of sigma2
  soma1.2 = matrix(0,p,qq)                  # Cross-derived: beta_ff
  soma1.3 = matrix(0,p,1)                   # Cross-derived: beta_sigma2
  soma2.3 = matrix(0,qq,1)                  # Cross-derived: ff_sigma2

  DeltaWei     = matrix(0,p+qq+1,n)         # For measures scheme case weigth
  DeltaSca     = matrix(0,p+qq+1,n)         # For measures scheme case scale
  DeltaRes     = matrix(0,p+qq+1,n)         # For measures scheme case resposta
  DeltaExp     = matrix(0,p+qq+1,n)         # For measures scheme case explanatory
  Score_vector = matrix(0,p+qq+1,n)         # For Score vector in case deletion measures

  for(i in 1:n){

    xi       = matrix(x[i,  ],ncol=p)
    Ni       = matrix(N[i,  ],ncol=qq)
    Deltai   = cbind(xi,Ni)
    mui      = Deltai%*%alpha

    soma1.1 = soma1.1 + u0[i]*(t(xi)%*%xi)                  # for beta beta
    soma2.2 = soma2.2 + u0[i]*(t(Ni)%*%Ni)                  # for ff ff
    soma3.3 = soma3.3 + (u2[i]-2*u1[i]*mui+u0[i]*mui^2)     # for sigma2 sigma2
    soma1.2 = soma1.2 + u0[i]*(t(xi)%*%Ni)                  # for beta ff
    soma1.3 = soma1.3 + t(xi)*as.numeric(u1[i]-mui*u0[i])   # for beta sigma2
    soma2.3 = soma2.3 + t(Ni)*as.numeric(u1[i]-mui*u0[i])   # for ff sigma2


    #======= Derivates for weight perturbation

    DeltaWbeta   = (1/sigma2) * (t(xi)%*%(u1[i]-u0[i]*mui))
    DeltaWff     = (1/sigma2) * (t(Ni)%*%(u1[i]-u0[i]*mui))
    DeltaWsigma2 = (-0.5/(sigma2)) + (0.5/(sigma2^2))*(u2[i]-2*u1[i]*mui+mui^2*u0[i])


    #======= Measures for case weight perturbation


    DeltaWei[,i] = c(t(DeltaWbeta), t(DeltaWff), t(DeltaWsigma2))


    #======= Derivates for scale perturbation

    DeltaSbeta   = (1/sigma2) * (t(xi)%*%(u1[i]-u0[i]*mui))
    DeltaSff     = (1/sigma2) * (t(Ni)%*%(u1[i]-u0[i]*mui))
    DeltaSsigma2 = (0.5/(sigma2^2))*(u2[i]-2*u1[i]*mui+mui^2*u0[i])


    #======= Measures for case scale


    DeltaSca[,i] = c(t(DeltaSbeta), t(DeltaSff), t(DeltaSsigma2))


    #======= Derivates for response variable perturbation

    DeltaRbeta   = -(sd(y)/sigma2) *  u0[i]*t(xi)
    DeltaRff     = -(sd(y)/sigma2) *  u0[i]*t(Ni)
    DeltaRsigma2 = (sd(y)/(sigma2^2)) * (u0[i]*mui-u1[i])

    #======= Measures for case response variable

    DeltaRes[,i] = c(t(DeltaRbeta), t(DeltaRff), t(DeltaRsigma2))


    #======= Derivates for explanatory variable perturbation

    if(k!=0){

      I_Sr     = matrix(0,p,1)               # For standard error in k-th explanatory variable perturbation
      I_Sr[k,] = 1
      Sr       = sd(x[,k])

      DeltaEbeta   = (Sr/sigma2)*(u1[i]*I_Sr-u0[i]*(I_Sr%*%mui+t(xi)%*%t(I_Sr)%*%beta))
      DeltaEff     = (-Sr/sigma2)*(u0[i]*t(Ni)%*%t(I_Sr)%*%beta)
      DeltaEsigma2 = (Sr/sigma2^2) * (u0[i]*(t(I_Sr)%*%beta%*%mui) - u1[i]* (t(I_Sr)%*%beta))

    }

    #======= Measures for case explanatory variable

    DeltaExp[,i] = c(t(DeltaEbeta), t(DeltaEff), t(DeltaEsigma2))


    #======= Derivates for case deletion measures: Score vector

    d.beta   = (1/sigma2)*(t(xi)%*%(u1[i]-u0[i]*mui))
    d.ff     = (1/sigma2)*(t(Ni)%*%(u1[i]-u0[i]*mui)) - (lambda/n)*K%*%ff
    d.sigma2 = (-0.5/sigma2) + (0.5/sigma2^2)*(u2[i]-2*u1[i]*mui+mui^2*u0[i])

    Score_vector[,i] = c(t(d.beta),t(d.ff), d.sigma2)

  }

  #==== Derivatives with respect to beta, ff and sigma2

  d.betabeta     = ((-1/sigma2)*soma1.1)
  d.ff           = ((-1/(sigma2))*soma2.2)-lambda*K
  d.sigma2sigma2 = (n/(2*sigma2^2))-((1/(sigma2^3))* soma3.3)
  d.betaff       = ((-1/(sigma2))*soma1.2)
  d.betasigma2   = ((-1/(sigma2^2))*soma1.3)
  d.fsigma2      = ((-1/(sigma2^2))*soma2.3)

  D.beta         = cbind(d.betabeta,d.betaff,d.betasigma2)
  D.f            = cbind(t(d.betaff),d.ff,d.fsigma2)
  D.sigma2       = cbind(t(d.betasigma2),t(d.fsigma2),d.sigma2sigma2)

  #=== Hessian matrix

  Q = cbind(t(D.beta), t(D.f), t(D.sigma2))

  #=== Inverse of the -Hessian matrix

  Inv_Q = solve(-Q)

  #=== Approximation theta [-i] for case deletion measures

  GD=matrix(0,nrow = n)
  GD1=matrix(0,nrow = n)

  for(j in 1:n){
    Score_vectorj = matrix(Score_vector[,j],nrow=(p+qq+1),ncol=1)
    theta_i       = theta + (Inv_Q%*%Score_vectorj)

    GD[j]  = t(theta_i-theta)%*%(-Q)%*%(theta_i-theta)
    GD1[j] = t(Score_vectorj)%*%Inv_Q%*%Score_vectorj

  }

  #===============================================================================

  #=== Detecting inluential observations, case weight -Q_w0

  Inf1 = t(DeltaWei)%*%Inv_Q%*%DeltaWei

  #=== Conformal normal curvature:

  Mo1  = diag(Inf1) /(sum(diag(Inf1)))

  #===============================================================================

  #=== Detecting inluential observations, case scale -Q_w0

  Inf2 = t(DeltaSca)%*%Inv_Q%*%DeltaSca

  #=== Conformal normal curvature:

  Mo2  = diag(Inf2) /(sum(diag(Inf2)))

  #===============================================================================

  #=== Detecting inluential observations, case response variable  -Q_w0

  Inf3 = t(DeltaRes)%*%Inv_Q%*%DeltaRes

  #=== Conformal normal curvature:

  Mo3  =  diag(Inf3) /(sum(diag(Inf3)))

  #===============================================================================

  #=== Detecting inluential observations, case explanatory variable  -Q_w0

  Inf4 = t(DeltaExp)%*%Inv_Q%*%DeltaExp

  #=== Conformal normal curvature:

  Mo4  =  diag(Inf4) /(sum(diag(Inf4)))

  #===============================================================================

    return(list(Hessian=Q, GD1=GD1, Curvature_W=Mo1,Curvature_S=Mo2,Curvature_R=Mo3,Curvature_E=Mo4))

}






