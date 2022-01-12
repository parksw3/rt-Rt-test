renewal_det <- function(Rfun=function(t) 2.5,
                        r1, r2,
                        theta=1.5,
                        N=40000,
                        dt=0.025,
                        genfun1=function(x) dgamma(x, 5, 1),
                        genfun2=function(x) dgamma(x, 5, 1),
                        genfun3=function(x) dgamma(x, 5, 1),
                        I01=0.1,
                        I02=0.001,
                        tmax=100,
                        genmax=2000) {
  gen1 <- genfun1(0:genmax*dt+dt)
  gen1 <- gen1/sum(gen1)
  gen2 <- genfun2(0:genmax*dt+dt)
  gen2 <- gen2/sum(gen2)
  gen3 <- genfun3(0:genmax*dt+dt)
  gen3 <- gen3/sum(gen3)
  tvec <- seq(0, tmax, by=dt)
  tveclong <- seq(-genmax*dt, tmax, by=dt)
  Ivec1 <- rep(0, tmax/dt+genmax)
  Ivec2 <- rep(0, tmax/dt+genmax)
  
  if (missing(r1)) {
    r1 <- optim(0, function(x) (1/sum(exp(-x*(0:genmax*dt+dt)) * gen1)-Rfun(0))^2,
                lower=-0.4,
                upper=0.4,
                method="Brent")[[1]]
    
  }
  
  if (missing(r2)) {
    r2 <- optim(0, function(x) (1/sum(exp(-x*(0:genmax*dt+dt)) * gen2)-theta*Rfun(0))^2,
                lower=-0.4,
                upper=0.4,
                method="Brent")[[1]]
    
  }
  
  Ivec1[1:(genmax+1)] <- exp(r1*1:(genmax+1)*dt)
  Ivec1[1:(genmax+1)] <- I01*Ivec1[1:(genmax+1)]/Ivec1[(genmax+1)]*dt
  
  Ivec2[1:(genmax+1)] <- exp(r2*1:(genmax+1)*dt)
  Ivec2[1:(genmax+1)] <- I02*Ivec2[1:(genmax+1)]/Ivec2[(genmax+1)]*dt
  
  Svec <- rep(0, tmax/dt)
  Svec[1] <- N - sum(Ivec1[1:(genmax+1)]) - sum(Ivec2[1:(genmax+1)])
  
  for (i in 2:length(tvec)) {
    j <- i + genmax
    
    Ivec1[j] <- Svec[i-1] * Rfun(tvec[i]) * sum(Ivec1[max(1, j-genmax):(j-1)] * gen1[(genmax+1):2])/N
    Ivec2[j] <- Svec[i-1] * theta * Rfun(tvec[i]) * sum(Ivec2[max(1, j-genmax):(j-1)] * gen2[(genmax+1):2])/N
    Svec[i] <- Svec[i-1] - Ivec1[i] - Ivec2[i]
  }
  
  Rt1 <- tail(sapply(tvec, Rfun), -1) * head(Svec, -1)/N
  Rt2 <- tail(sapply(tvec, Rfun), -1) * head(Svec, -1)/N * theta
  
  Rtest1 <- tail(Ivec1, -1)/sapply(2:(length(Ivec1)), function(x) sum(Ivec1[max(1, x-genmax):(x-1)] * gen3[min(x, genmax+1):2]))
  Rtest2 <- tail(Ivec2, -1)/sapply(2:(length(Ivec2)), function(x) sum(Ivec2[max(1, x-genmax):(x-1)] * gen3[min(x, genmax+1):2]))
  
  Rtest1 <- Rtest1[-c(1:genmax)]
  Rtest2 <- Rtest2[-c(1:genmax)]
  
  Ivec1 <- Ivec1[-c(1:genmax)]
  Ivec2 <- Ivec2[-c(1:genmax)]
  
  data.frame(
    tvec=tvec[-1],
    Ivec1=Ivec1[-1]/dt,
    Ivec2=Ivec2[-1]/dt,
    Svec=Svec[-1],
    Rt1=Rt1,
    Rt2=Rt2,
    Rtest1=Rtest1,
    Rtest2=Rtest2,
    r1=c(diff(Ivec1))/head(Ivec1, -1)/dt,
    r2=c(diff(Ivec2))/head(Ivec2, -1)/dt
  )
}
