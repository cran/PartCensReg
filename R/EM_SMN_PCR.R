################################################################################
## Observed information matrix / EM-type algorithm  for the SMN-PCR model ##
################################################################################

funcIi <- function(w, nu, di, type){
  if(type=="Normal"){
    I_phi   <- sqrt(2*pi)*dnorm(sqrt(di),0,1)
    I_Acphi <- (1/sqrt(2*pi))*I_phi
  }
  if(type=="T"){
    I_phi   <- (nu/2)^((nu/2))*(gamma((nu/2)+w)/gamma(nu/2))*((nu+di)/2)^(-((nu/2)+w))
    I_Acphi <- (1/sqrt(2*pi))*I_phi
  }
  if(type=="Slash"){
    func    <- function(u) nu*u^(nu+w-1)*exp(-0.5%*%di%*%u)
    resp    <- integrate(func,0,1)$value
    I_phi   <- resp
    I_Acphi <- (1/sqrt(2*pi))*I_phi
  }
  if(type=="NormalC"){
    I_phi   <- sqrt(2*pi)*(nu[2]*(nu[1]^(w-0.5))*dnorm(sqrt(di),0,sqrt(1/nu[1]))) + sqrt(2*pi)*((1-nu[2])*dnorm(sqrt(di),0,1))
    I_Acphi <- (1/sqrt(2*pi))*I_phi
  }
  I_phi   <- as.numeric(I_phi)
  I_Acphi <- as.numeric(I_Acphi)
  return(list(I_phi=I_phi, I_Acphi=I_Acphi))
}

funcKi <- function(yi, mui, sigma2, nu, type){
  aux  <- (yi-mui)/sqrt(sigma2)
  if(type=="Normal") {Ki <- sqrt(2*pi*sigma2)*dnorm(aux,0,1)/sqrt(sigma2)}
  if(type=="T")      {Ki <- sqrt(2*pi*sigma2)*dt(aux,df=nu)/sqrt(sigma2)}
  if(type=="Slash")  {Ki <- sqrt(2*pi*sigma2)*dSlash(aux,0,1,nu)/sqrt(sigma2)}
  if(type=="NormalC"){Ki <- sqrt(2*pi*sigma2)*dNormalC(aux,0,1,nu)/sqrt(sigma2)}
  Ki <- as.numeric(Ki)
  return(list(Ki=Ki))
}

funcJi <- function(yi, mui, sigma2, nu, type){
  aux <- (yi-mui)/sqrt(sigma2)
  if(type=="Normal") {Ji <- pnorm(aux)}
  if(type=="T")      {Ji <- pt(aux,df=nu)}
  if(type=="Slash")  {Ji <- AcumSlash(aux,0,1,nu)}
  if(type=="NormalC"){Ji <- AcumNormalC(aux,0,1,nu)}
  Ji <- as.numeric(Ji)
  return(list(Ji=Ji))
}


MIobsSMN <- function(x, y, c, cens, tt, nu, beta, ff, sigma2, lambda, type){

  p  <- ncol(x)
  n  <- nrow(x)

  tr <- sort(unique(tt))
  r  <- length(tr)
  t2 <- ncs(tr, nknots=r)
  K  <- attr(t2, "K")
  N  <- matrix(data = 0,nrow = n,ncol = r)
  for(i in 1:r){N[tt==tr[i],i] = 1}
  qq <-ncol(N)
  Ks <- bdiag(matrix(0,p,p),K)

  alpha <- as.vector(c(beta,ff))

  if(cens=="1"){

    soma1.2 <- matrix(0,p+qq,p+qq)
    soma2.2 <- matrix(0,p+qq,1)
    soma3.2 <- 0
    soma1.3 <- matrix(0,p+qq,p+qq)
    soma2.3 <- matrix(0,p+qq,1)
    soma3.3 <- 0

    for(i in 1:n){

      yi     <- y[i]
      xi     <- matrix(x[i,  ],ncol=p)
      Ni     <- matrix(N[i,  ],ncol=qq)
      Deltai <- cbind(xi,Ni)
      mui    <- Deltai%*%alpha

      if(c[i]==0){

        Kti         <- funcKi(yi, mui, sigma2, nu, type)$Ki
        di          <- (yi-mui)^2/sigma2
        Iphi_3.2    <- funcIi(3/2, nu, di, type)$I_phi
        Iphi_5.2    <- funcIi(5/2, nu, di, type)$I_phi

        di_gama     <- (-2/sigma2)*(t(Deltai)%*%(yi-mui))
        di_sig      <- as.numeric((-1/sigma2^2)*(yi-mui)^2)
        di_gamagama <- (2/sigma2)*(t(Deltai)%*%Deltai)
        di_sigsig   <- as.numeric((2/sigma2^3)*(yi-mui)^2)
        di_gamasig  <- (2/sigma2^2)*(t(Deltai)%*%(yi-mui))

        dK_gama     <- (-0.5)*Iphi_3.2*di_gama
        dK_sig      <-  (-0.5)*Iphi_3.2*di_sig
        dK_gamagama <- (-0.5)*Iphi_3.2*di_gamagama + 0.25*Iphi_5.2*di_gama%*%t(di_gama)
        dK_sigsig   <- (-0.5)*Iphi_3.2*di_sigsig + 0.25*Iphi_5.2*di_sig*di_sig
        dK_gamasig  <- (-0.5)*Iphi_3.2*di_gamasig + 0.25*Iphi_5.2*di_sig*di_gama

        soma1.2 <- soma1.2 + (((1/Kti^2)*dK_gama%*%t(dK_gama))-((1/Kti)*dK_gamagama))
        soma2.2 <- soma2.2 + (((1/Kti^2)*dK_sig*dK_gama)-((1/Kti)*dK_gamasig))
        soma3.2 <- soma3.2 + (((1/Kti^2)*dK_sig*dK_sig)-((1/Kti)*dK_sigsig))
      }

      if(c[i]==1){

        Jti         <- funcJi(yi, mui, sigma2, nu, type)$Ji
        di          <- (yi-mui)^2/sigma2
        Ai          <- (yi-mui)/sqrt(sigma2)
        IAcphi_1.2  <- funcIi(1/2, nu, di, type)$I_Acphi
        IAcphi_3.2  <- funcIi(3/2, nu, di, type)$I_Acphi

        di_gama     <- (-2/sigma2)*(t(Deltai)%*%(yi-mui))
        di_sig      <- as.numeric((-1/sigma2^2)*(yi-mui)^2)
        Ai_gama     <- (-1/sqrt(sigma2))*t(Deltai)
        Ai_sig      <- as.numeric((-0.5/sigma2^(3/2))*(yi-mui))
        Ai_gamagama <- matrix(0,p+qq,p+qq)
        Ai_sigsig   <- as.numeric((0.75/sigma2^(5/2))*(yi-mui))
        Ai_gamasig  <- (0.5/sigma2^(3/2))*t(Deltai)

        dJ_gama     <- IAcphi_1.2*Ai_gama
        dJ_sig      <-  IAcphi_1.2*Ai_sig
        dJ_gamagama <- (-0.5)*IAcphi_3.2*di_gama%*%t(Ai_gama) + IAcphi_1.2*Ai_gamagama
        dJ_sigsig   <- (-0.5)*IAcphi_3.2*di_sig*Ai_sig + IAcphi_1.2*Ai_sigsig
        dJ_gamasig  <- (-0.5)*IAcphi_3.2*di_sig*Ai_gama + IAcphi_1.2*Ai_gamasig

        soma1.3 <- soma1.3 + (((1/Jti^2)*dJ_gama%*%t(dJ_gama))-((1/Jti)*dJ_gamagama))
        soma2.3 <- soma2.3 + (((1/Jti^2)*dJ_sig*dJ_gama)-((1/Jti)*dJ_gamasig))
        soma3.3 <- soma3.3 + (((1/Jti^2)*dJ_sig*dJ_sig)-((1/Jti)*dJ_sigsig))
      }
    }

    I1_gamagama <- lambda*Ks
    I1_gamasig  <- matrix(0,p+qq,1)
    I1_sigsig   <- (-0.5/sigma2^2)*(n-sum(c))
    I2_gamagama <- soma1.2
    I2_gamasig  <- soma2.2
    I2_sigsig   <- soma3.2
    I3_gamagama <- soma1.3
    I3_gamasig  <- soma2.3
    I3_sigsig   <- soma3.3

    I1 <- rbind(cbind(I1_gamagama,I1_gamasig),cbind(t(I1_gamasig),I1_sigsig))
    I2 <- rbind(cbind(I2_gamagama,I2_gamasig),cbind(t(I2_gamasig),I2_sigsig))
    I3 <- rbind(cbind(I3_gamagama,I3_gamasig),cbind(t(I3_gamasig),I3_sigsig))

    MI <- I1 + I2 + I3

  }

  if(cens=="2"){

    soma1.2 <- matrix(0,p+qq,p+qq)
    soma2.2 <- matrix(0,p+qq,1)
    soma3.2 <- 0
    soma1.3 <- matrix(0,p+qq,p+qq)
    soma2.3 <- matrix(0,p+qq,1)
    soma3.3 <- 0

    for(i in 1:n){

      yi    <- y[i]
      xi    <- matrix(x[i,  ],ncol=p)
      Ni    <- matrix(N[i,  ],ncol=qq)
      Deltai<- cbind(xi,Ni)
      mui   <- Deltai%*%alpha

      if(c[i]==0){

        Kti         <- funcKi(yi, mui, sigma2, nu, type)$Ki
        di          <- (yi-mui)^2/sigma2
        Iphi_3.2    <- funcIi(3/2, nu, di, type)$I_phi
        Iphi_5.2    <- funcIi(5/2, nu, di, type)$I_phi

        di_gama     <- (-2/sigma2)*(t(Deltai)%*%(yi-mui))
        di_sig      <- as.numeric((-1/sigma2^2)*(yi-mui)^2)
        di_gamagama <- (2/sigma2)*(t(Deltai)%*%Deltai)
        di_sigsig   <- as.numeric((2/sigma2^3)*(yi-mui)^2)
        di_gamasig  <- (2/sigma2^2)*(t(Deltai)%*%(yi-mui))

        dK_gama     <- (-0.5)*Iphi_3.2*di_gama
        dK_sig      <-  (-0.5)*Iphi_3.2*di_sig
        dK_gamagama <- (-0.5)*Iphi_3.2*di_gamagama + 0.25*Iphi_5.2*di_gama%*%t(di_gama)
        dK_sigsig   <- (-0.5)*Iphi_3.2*di_sigsig + 0.25*Iphi_5.2*di_sig*di_sig
        dK_gamasig  <- (-0.5)*Iphi_3.2*di_gamasig + 0.25*Iphi_5.2*di_sig*di_gama

        soma1.2 <- soma1.2 + (((1/Kti^2)*dK_gama%*%t(dK_gama))-((1/Kti)*dK_gamagama))
        soma2.2 <- soma2.2 + (((1/Kti^2)*dK_sig*dK_gama)-((1/Kti)*dK_gamasig))
        soma3.2 <- soma3.2 + (((1/Kti^2)*dK_sig*dK_sig)-((1/Kti)*dK_sigsig))
      }

      if(c[i]==1){

        Jti         <- funcJi(yi, mui, sigma2, nu, type)$Ji
        di          <- (yi-mui)^2/sigma2
        Ai          <- (yi-mui)/sqrt(sigma2)
        IAcphi_1.2  <- funcIi(1/2, nu, di, type)$I_Acphi
        IAcphi_3.2  <- funcIi(3/2, nu, di, type)$I_Acphi

        di_gama     <- (-2/sigma2)*(t(Deltai)%*%(yi-mui))
        di_sig      <- as.numeric((-1/sigma2^2)*(yi-mui)^2)
        Ai_gama     <- (-1/sqrt(sigma2))*t(Deltai)
        Ai_sig      <- as.numeric((-0.5/sigma2^(3/2))*(yi-mui))
        Ai_gamagama <- matrix(0,p+qq,p+qq)
        Ai_sigsig   <- as.numeric((0.75/sigma2^(5/2))*(yi-mui))
        Ai_gamasig  <- (0.5/sigma2^(3/2))*t(Deltai)

        dJ_gama     <- IAcphi_1.2*Ai_gama
        dJ_sig      <-  IAcphi_1.2*Ai_sig
        dJ_gamagama <- (-0.5)*IAcphi_3.2*di_gama%*%t(Ai_gama) + IAcphi_1.2*Ai_gamagama
        dJ_sigsig   <- (-0.5)*IAcphi_3.2*di_sig*Ai_sig + IAcphi_1.2*Ai_sigsig
        dJ_gamasig  <- (-0.5)*IAcphi_3.2*di_sig*Ai_gama + IAcphi_1.2*Ai_gamasig

        soma1.3 <- soma1.3 + (((1/(1-Jti)^2)*dJ_gama%*%t(dJ_gama))+((1/(1-Jti))*dJ_gamagama))
        soma2.3 <- soma2.3 + (((1/(1-Jti)^2)*dJ_sig*dJ_gama)+((1/(1-Jti))*dJ_gamasig))
        soma3.3 <- soma3.3 + (((1/(1-Jti)^2)*dJ_sig*dJ_sig)+((1/(1-Jti))*dJ_sigsig))
      }
    }

    I1_gamagama <- lambda*Ks
    I1_gamasig  <- matrix(0,p+qq,1)
    I1_sigsig   <- (-0.5/sigma2^2)*(n-sum(c))
    I2_gamagama <- soma1.2
    I2_gamasig  <- soma2.2
    I2_sigsig   <- soma3.2
    I3_gamagama <- soma1.3
    I3_gamasig  <- soma2.3
    I3_sigsig   <- soma3.3

    I1 <- rbind(cbind(I1_gamagama,I1_gamasig),cbind(t(I1_gamasig),I1_sigsig))
    I2 <- rbind(cbind(I2_gamagama,I2_gamasig),cbind(t(I2_gamasig),I2_sigsig))
    I3 <- rbind(cbind(I3_gamagama,I3_gamasig),cbind(t(I3_gamasig),I3_sigsig))

    MI <- I1 + I2 + I3

  }

  return(list(MI=MI))

}

EMSpline.censFinal_MobsSMN <- function(x, y, c, cens, tt, nu, error, iter.max, type, delta.in=NA, lambda.FIX, nu.FIX, lambda.in, k=0){

  if(nu.FIX==FALSE){

  p   <- ncol(x)
  n   <- nrow(x)
  reg <- lm(y ~ x -1)

  tr  <- sort(unique(tt))
  r   <- length(tr)
  t2  <- ncs(tr, nknots=r)
  K   <- attr(t2, "K")
  N   <- matrix(data = 0,nrow = n,ncol = r)
  for(i in 1:r){N[tt==tr[i],i] = 1}
  qq  <- ncol(N)

  lambda <- lambda.in

  #=== Initial values

  beta   <- as.vector(coefficients(reg),mode="numeric")
  sigma2 <- as.numeric(sum((y-(x%*%beta))^2)/(n-p))
  ff     <- solve(t(N)%*%N+lambda*sigma2*K)%*%t(N)%*%(y-x%*%beta)
  Delta  <- cbind(x,N)
  alpha  <- as.vector(c(beta,ff))
  mu     <- Delta%*%alpha
  delta  <- delta.in

  theta  <- c(beta,sigma2,nu,lambda)

  if(cens=="1"){Lim1 <- rep(-Inf,n); Lim2 <- y}
  if(cens=="2"){Lim1 <- y; Lim2 <- rep(Inf,n)}

  if(type=="Normal") {ver_velho <- sum(log(dNormal(c, y, mu, sigma2, cens)))}
  if(type=="T")      {ver_velho <- sum(log(dT(c, y, mu, sigma2, nu, cens)))}
  if(type=="Slash")  {ver_velho <- sum(log(dSL(c, y, mu, sigma2, nu, cens)))}
  if(type=="NormalC"){ver_velho <- sum(log(dCN(c, y, mu, sigma2, nu, cens)))}

  criterio <- 1
  count    <- 0

  while((criterio > error) & (count <= iter.max)){

    pb <- txtProgressBar(min = 1,max = iter.max,style = 3)

    NCensEUY     <- NCensurEsperUY(y,mu,sigma2,nu,delta,type)
    u0           <- NCensEUY$EUY0
    u1           <- NCensEUY$EUY1
    u2           <- NCensEUY$EUY2

    if(sum(c)>0){
      CensEUY    <- CensEsperUY1(mu[c==1],sigma2=sigma2,nu=nu,delta,Lim1=Lim1[c==1],Lim2=Lim2[c==1],type,cens)
      u0[c==1]   <- CensEUY$EUY0
      u1[c==1]   <- CensEUY$EUY1
      u2[c==1]   <- CensEUY$EUY2
    }

    soma1 <- matrix(0,p,p)
    soma2 <- matrix(0,p,1)
    soma3 <- matrix(0,qq,qq)
    soma4 <- matrix(0,qq,1)
    soma5 <- 0

    for(i in 1: n){
      xi     <- matrix(x[i,  ],ncol=p)
      Ni     <- matrix(N[i,  ],ncol=qq)
      Deltai <- cbind(xi,Ni)
      mui    <- Deltai%*%alpha
      soma1  <- soma1 + u0[i]*t(xi)%*%xi
      soma2  <- soma2 + t(xi)%*%(u1[i]-u0[i]*(Ni%*%ff))
      soma3  <- soma3 + u0[i]*t(Ni)%*%Ni
      soma4  <- soma4 + t(Ni)%*%(u1[i]-u0[i]*(xi%*%beta))
      soma5  <- soma5 + u2[i]-2*u1[i]*mui+mui^2*u0[i]
    }

    beta   <- solve(soma1)%*%soma2
    ff     <- solve(soma3+(lambda*sigma2*K))%*%soma4
    alpha  <- rbind(beta,ff)
    sigma2 <- as.numeric(soma5/n)

    mu     <- Delta%*%alpha

    if(type=="Normal"){ver_novo <- sum(log(dNormal(c, y, mu, sigma2, cens)))}
    if(type=="T"){
      f  <- function(nu){sum(log(dT(c, y, mu, sigma2, nu, cens)))}
      nu <- optimize(f, c(2.01,150), tol = 0.000001, maximum = TRUE)$maximum
      ver_novo <- f(nu)}
    if(type=="Slash"){
      f  <- function(nu){sum(log(dSL(c, y, mu, sigma2, nu, cens)))}
      nu <- optimize(f, c(0.01,15), tol = 0.000001, maximum = TRUE)$maximum
      ver_novo <- f(nu)}
    if(type=="NormalC"){
      nu1<-nu[1];nu2<-nu[2];
      f<-function(nu){
        nu1<-nu[1]
        nu2<-nu[2]
        return(-sum(log(dCN(c, y, mu, sigma2,c(nu1,nu2), cens))))
      }
      nus <- optimx(c(nu1,nu2), method = "nlminb", fn=f, lower=c(0.1,0.1), upper=c(0.9,0.9), hessian = TRUE,control=list(maximize=FALSE))
      nu  <- c(nus$p1,nus$p2)
      ver_novo <- -f(nu)}

    ft <- function(lambda) AICsmn(lambda,ver_novo,ff,K,N,sigma2,p,type)
    if(lambda.FIX==FALSE){lambda<- optimize(f=ft, interval=c(0,100),maximum=TRUE,tol=10^{-6})$maximum}

    theta     <- c(beta,sigma2,nu,lambda)
    criterio  <- (abs(1-ver_novo/ver_velho))
    ver_velho <- ver_novo

    count <- (count+1)
    AIC<- -ft(lambda)
    setTxtProgressBar(pb, count)
  }

  MIf <- MIobsSMN(x, y, c, cens, tt, nu, beta, ff, sigma2, lambda, type)
  MI  <- MIf$MI
  Diagnostics <- dignostic.EM_Sem(x, y, tt,  beta, ff, sigma2, lambda,  u0, u1, u2, k, type)
  close(pb)
  return(list(beta=beta,sigma2=sigma2,lambda=lambda, AIC= AIC, ff=ff, yest= mu, loglik = ver_novo, iter=count, nu=nu, theta=theta, MI=MI, u0=u0, u1=u1, u2=u2, D=Diagnostics))

  }


  if(nu.FIX==TRUE){

    p   <- ncol(x)
    n   <- nrow(x)
    reg <- lm(y ~ x -1)

    tr  <- sort(unique(tt))
    r   <- length(tr)
    t2  <- ncs(tr, nknots=r)
    K   <- attr(t2, "K")
    N   <- matrix(data = 0,nrow = n,ncol = r)
    for(i in 1:r){N[tt==tr[i],i] = 1}
    qq  <- ncol(N)

    lambda <- lambda.in

    #=== Initial values

    beta   <- as.vector(coefficients(reg),mode="numeric")
    sigma2 <- as.numeric(sum((y-(x%*%beta))^2)/(n-p))
    ff     <- solve(t(N)%*%N+lambda*sigma2*K)%*%t(N)%*%(y-x%*%beta)
    Delta  <- cbind(x,N)
    alpha  <- as.vector(c(beta,ff))
    mu     <- Delta%*%alpha
    delta  <- delta.in

    theta  <- c(beta,sigma2,nu,lambda)

    if(cens=="1"){Lim1 <- rep(-Inf,n); Lim2 <- y}
    if(cens=="2"){Lim1 <- y; Lim2 <- rep(Inf,n)}

    if(type=="Normal") {ver_velho <- sum(log(dNormal(c, y, mu, sigma2, cens)))}
    if(type=="T")      {ver_velho <- sum(log(dT(c, y, mu, sigma2, nu, cens)))}
    if(type=="Slash")  {ver_velho <- sum(log(dSL(c, y, mu, sigma2, nu, cens)))}
    if(type=="NormalC"){ver_velho <- sum(log(dCN(c, y, mu, sigma2, nu, cens)))}

    criterio <- 1
    count    <- 0

    while((criterio > error) & (count <= iter.max)){

      pb <- txtProgressBar(min = 1,max = iter.max,style = 3)

      NCensEUY     <- NCensurEsperUY(y,mu,sigma2,nu,delta,type)
      u0           <- NCensEUY$EUY0
      u1           <- NCensEUY$EUY1
      u2           <- NCensEUY$EUY2

      if(sum(c)>0){
        CensEUY    <- CensEsperUY1(mu[c==1],sigma2=sigma2,nu=nu,delta,Lim1=Lim1[c==1],Lim2=Lim2[c==1],type,cens)
        u0[c==1]   <- CensEUY$EUY0
        u1[c==1]   <- CensEUY$EUY1
        u2[c==1]   <- CensEUY$EUY2
      }

      soma1 <- matrix(0,p,p)
      soma2 <- matrix(0,p,1)
      soma3 <- matrix(0,qq,qq)
      soma4 <- matrix(0,qq,1)
      soma5 <- 0

      for(i in 1: n){
        xi     <- matrix(x[i,  ],ncol=p)
        Ni     <- matrix(N[i,  ],ncol=qq)
        Deltai <- cbind(xi,Ni)
        mui    <- Deltai%*%alpha
        soma1  <- soma1 + u0[i]*t(xi)%*%xi
        soma2  <- soma2 + t(xi)%*%(u1[i]-u0[i]*(Ni%*%ff))
        soma3  <- soma3 + u0[i]*t(Ni)%*%Ni
        soma4  <- soma4 + t(Ni)%*%(u1[i]-u0[i]*(xi%*%beta))
        soma5  <- soma5 + u2[i]-2*u1[i]*mui+mui^2*u0[i]
      }

      beta   <- solve(soma1)%*%soma2
      ff     <- solve(soma3+(lambda*sigma2*K))%*%soma4
      alpha  <- rbind(beta,ff)
      sigma2 <- as.numeric(soma5/n)

      mu     <- Delta%*%alpha

      if(type=="Normal"){ver_novo <- sum(log(dNormal(c, y, mu, sigma2, cens)))}
      if(type=="T"){
        f  <- function(nu){sum(log(dT(c, y, mu, sigma2, nu, cens)))}
        nu=nu
        ver_novo <- f(nu)}
      if(type=="Slash"){
        f  <- function(nu){sum(log(dSL(c, y, mu, sigma2, nu, cens)))}
        nu=nu
        ver_novo <- f(nu)}
      if(type=="NormalC"){
        nu1<-nu[1];nu2<-nu[2];
        f<-function(nu){
          nu1<-nu[1]
          nu2<-nu[2]
          return(sum(log(dCN(c, y, mu, sigma2,c(nu1,nu2), cens))))
        }
        ver_novo <- f(nu)}

      ft <- function(lambda) AICsmn(lambda,ver_novo,ff,K,N,sigma2,p,type)
      if(lambda.FIX==FALSE){lambda<- optimize(f=ft, interval=c(0,100),maximum=TRUE,tol=10^{-6})$maximum}

      theta     <- c(beta,sigma2,nu,lambda)
      criterio  <- (abs(1-ver_novo/ver_velho))
      ver_velho <- ver_novo

      count <- (count+1)
      AIC<- -ft(lambda)
      setTxtProgressBar(pb, count)
      }

    MIf <- MIobsSMN(x, y, c, cens, tt, nu, beta, ff, sigma2, lambda, type)
    MI  <- MIf$MI
    Diagnostics <- dignostic.EM_Sem(x, y, tt,  beta, ff, sigma2, lambda,  u0, u1, u2, k, type)
    close(pb)
    return(list(beta=beta,sigma2=sigma2,lambda=lambda, AIC= AIC, ff=ff, yest= mu, loglik = ver_novo, iter=count, nu=nu, MI=MI,D=Diagnostics))

  }

}

