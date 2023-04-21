rgld.fit<-function (xdat, r = dim(xdat)[2],n=20,ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,  hl= NULL, 
                             mulink = identity, siglink = identity, shlink = identity, hlink = identity, show = TRUE, 
                             method = "Nelder-Mead", maxit = 10000, ...) 
{
  
  options(digits=8)
  z <- list()
  k <- list()
  
  npmu <- length(mul) + 1
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  
  if (is.null(mul)) {
    mumat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
  }
  
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
  }
  
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, dim(xdat)[1]))
  } else {
    shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
  }
  
  
  init <- ginit3(parglo(lmoms(xdat[,1]))$para,n=n)
  
  z$link <- deparse(substitute(c(mulink, siglink, shlink)))
  
  z1 <- apply(xdat, 1, max, na.rm = TRUE) # z_1
  zr <- apply(xdat, 1, min, na.rm = TRUE) # z_r
  
  
  gld.lik <- function(a) {
    
    mu <- mulink(mumat %*% (a[1:npmu]))
    sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y
    
    f <- 1 + y^(1/xi)
    
    F_ <- 1-(1-1/f)
    
    if (any(y <= 0,na.rm=T) || any(sc <= 0,na.rm=T) || any(f <=0,na.rm=T) || any(F_ >1,na.rm=T))  
      return(10^6)
    sum(log(sc)) + sum(log(1 + y^(1/xi)) * 2 ) + sum(log(y) * (1 - 1/xi)) 
    
  }
  
  
  rgld.lik <- function(a) {
    
    
    mu <- mulink(drop(mumat %*% (a[1:npmu])))
    sc <- siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
    xi <- shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
    
    
    ri  <- (r-seq(1:(r))) # r-i
    cr  <- (1+ri)    # c_r
    
    # constraints 1 #
    
    if (any(sc <= 0,na.rm=T) | any(cr < 0,na.rm=T) ) return(10^6) 
    
    y <- 1 - xi * (xdat - mu)/sc
    f <- 1 + (1 - xi * (zr - mu)/sc)^(1/xi) 
    
    # constraints 2,3,4 #
    
    cond <-(sc[1]/xi[1])+mu[1]
    
    
    if (any(y<= 0, na.rm = TRUE)  || any(f<= 0, na.rm = TRUE) ) return(10^6)
    
    y <- log(sc) + (1 - 1/xi) * log(y) - log(cr)  
    y <- rowSums(y, na.rm = TRUE)
    
    sum((r + 1) * log(f) + y)
    
  }
  
  
  
  if(r==1){
    
    tryCatch(
      
      for(i in 1:nrow(init)){       
        
        value <- try(solnp(init[i,], gld.lik, LB =c(-Inf,0,-Inf) ,UB =c(Inf,Inf, Inf)))
        
        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
        }
        
      }
    )
    
  }else{
    
    tryCatch(
      
      for(i in 1:nrow(init)){       
        
        
        value <- try(solnp(init[i,], rgld.lik, LB =c(-Inf,0,-Inf) ,UB =c(Inf,Inf, Inf)))
        
        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
        }
        
      }
    )
    
  }
  
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, function(x) x$value[which.min(x$value)])) #%>% filter(value!=10^6)
  optim_value  <-optim_value[optim_value$value!=10^6,]
  
  if(r==1){optim_grad   <-sapply(optim_value$num, function(x) sum(abs(grad(gld.lik,k[[x]]$par))))
  }else   {optim_grad   <-sapply(optim_value$num, function(x) sum(abs(grad(rgld.lik,k[[x]]$par))))}
  
  optim_value$grad <- optim_grad
  
  optim_table1 <-optim_value[order(optim_value$grad, optim_value$value),]
  optim_table2 <-optim_table1[(optim_table1$grad>-1 & optim_table1$grad<1),]
  optim_table3 <-optim_table2[order(optim_table2$value),]
  
  selc_num  <- optim_table3[1,"num"]
  
  x  <-k[[selc_num]]
  
  mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
  sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
  xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
  
  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$data <- xdat
  z$mle  <- x$par
  
  if(r==1){z$grad <- grad(gld.lik,x$par)
  }else   {z$grad <- grad(rgld.lik,x$par)}
  
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$vals <- cbind(mu, sc, xi)
  z$r <- r
  
  z$rl20  <- gld.rl_m(z$mle,z$cov,year=20)
  z$rl50  <- gld.rl_m(z$mle,z$cov,year=50)
  z$rl100 <- gld.rl_m(z$mle,z$cov,year=100)
  z$rl200 <- gld.rl_m(z$mle,z$cov,year=200)
  # 
  nr   <-nrow(xdat) 
  
  options(digits=8)
  
  z$rslt <- round(data.frame(r = r,
                             nllh = z$nllh,
                             tr_V = tr(z$cov),
                             logdet = log(det(z$cov)),
                             aic  = cal_stat(z$nllh,3,method="AIC"),
                             bic  = cal_stat(z$nllh,3,nr,method="BIC"),
                             mle  = t(z$mle),
                             se   = t(z$se),
                             rl20   = z$rl20$rl,
                             rl20se = z$rl20$rl_se,
                             rl50   = z$rl50$rl,
                             rl50se = z$rl50$rl_se,
                             rl100   = z$rl100$rl,
                             rl100se = z$rl100$rl_se,
                             rl200   = z$rl200$rl,
                             rl200se = z$rl200$rl_se),3)
  
  
  z$g   <-  data.frame(r=r,
                       nllh=z$nllh,
                       g1 = z$grad[1],
                       g2 = z$grad[2],
                       g3 = z$grad[3]) 
  
  
  if(abs(sum(z$g[3:5]))>1) stop("grad is grater than 1")
  
  return(z)
  invisible(z)
}


gld.rl_m<-function(par, varcov, year){
  
  options(digits=22)
  z   <-list()
  
  del <-matrix(ncol=1,nrow=3)
  
  mu = par[1]
  sig= par[2]
  xi = par[3]
  
  f = 1-(1/year) # F
  f_inv = 1/f
  
  z$rl = mu +sig/xi - sig/xi * exp(xi*log(f_inv-1))
  
  del[1,1] <- 1
  
  del[2,1] <- (1-exp(xi*log(f_inv-1)))/xi
  
  del[3,1] <- -(sig/xi^2) * ( - exp(xi*log(f_inv-1)) * (log(f_inv-1)*xi-1) -1)
  
  del.t <-t(del)
  
  z$rl_se <-sqrt((del.t %*% varcov %*% del))
  
  invisible(z)
  
  
}


cal_stat<-function(x,npar,nr,method="AIC"){
  if(method =="AIC"){
    k=2
    value = 2*x+k*npar
  }else if(method =="BIC"){
    value = 2*x+log(nr)*npar
  }
  return(value)
}


ginit3<-function(par,n){
  
  init <-matrix(rep(0,n*3),ncol=3)
  
  init[1:n,1] <-par[1]
  init[1:n,2] <-par[2]
  init[1:n,3] <-par[3]
  
  
  init[2:n,1] <-par[1]+rnorm(n=n-1,mean=0,sd = 0.1)
  init[2:n,2] <-par[2]+rnorm(n=n-1,mean=0,sd = 0.1)
  init[2:n,3] <-par[3]+rnorm(n=n-1,mean=0,sd = 0.1)
  
  
  return(init)
  
}


qgld<-function (p,par) 
{
  
  mu = par[1]
  sig= par[2]
  xi = par[3]
  
  f = p
  f_inv = 1/f
  
  z = mu +sig/xi - sig/xi * exp(xi*log(f_inv-1))
  z
}
