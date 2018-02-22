## functions for fitting rr lyrae light curves

## fits RR Lyrae template
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##       omegas : vector of frequencies
##          tem : input templates, see make_template.R for a description
##           NN : number of newton steps at each frequency
##   use.errors : should photometric errors be used
##     use.dust : should dust (E[B-V]) be fit
##
##
## value
##          rss : the residual sum of squares at each frequency in omegas     
FitTemplate <- function(lc,omegas,tem,NN=5,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    ## center times at 0, makes phi a smoother function of omega
    lc[,1] <- lc[,1] - mean(lc[,1])
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    coeffs <- c(0,0,0,runif(1))
    ##rss_max <- sum((lm(m~dust,weights=weights)$residuals^2)*weights)
    rss <- rep(0,length(omegas))
    betas <- tem$abs_mag(1/omegas,tem) ## obtain absolute magnitudes at all frequencies
    for(ii in 1:length(omegas)){
        m_temp <- m - rep.int(betas[ii,],nb) ## correct for absolute magnitude
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omegas[ii],m_temp,t,dust,weights,nb,
                                   tem$template_funcs,tem$templated_funcs,use.errors,use.dust)
        }
        gammaf <- ConstructGamma(t,nb,coeffs[4],omegas[ii],tem$template_funcs)
        resid <- m_temp - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
        ##rss[ii] <- min(sum(weights*resid^2),rss_max)
        rss[ii] <- sum(weights*resid^2)
    }
    return(rss)
}

## for a given omega, return coefficients
## example usage: run FitTemplate to determine rss as function of omegas,
## find omega which minimizes rss, then find coeffs for this omega using ComputeCoeffs
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##        omega : frequency
##          tem : input templates
##           NN : number of newton steps, probably 10+ since no warm start
##   use.errors : should photometric errors be used
##     use.dust : should dust (E[B-V]) be fit
##
##
## value
##       coeffs : vector of [distance mod,ebv,peak-to-peak g amp,phase]
ComputeCoeffs <- function(lc,omega,tem,NN=20,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    ## center times at 0, makes phi a smoother function of omega
    mean_time <- mean(lc[,1])
    lc[,1] <- lc[,1] - mean_time
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    ##lc$mag <- lc$mag - rep.int(tem$betas,nb)
    m <- m - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    coeffs <- c(0,0,0,runif(1))
    J <- 0
    ## J prevents infinite loops
    while(coeffs[3]==0 & J < 10){
        for(jj in 1:NN){
            coeffs <- NewtonUpdate(coeffs[4],omega,m,t,dust,weights,nb,
                                   tem$template_funcs,tem$templated_funcs,use.errors,use.dust)
        }
        J <- J + 1
    }
    ## shift phase back to original scale
    coeffs[4] <- (coeffs[4] - omega*mean_time) %% 1
    return(coeffs)
}


## for a given omega, phi return coefficients mu, dust, amp
## example usage: run FitTemplate to determine rss as function of omegas,
## find omega which minimizes rss, then find coeffs for this omega using ComputeCoeffs
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##        omega : frequency
##          phi : phase
##          tem : input templates
##   use.errors : should photometric errors be used
##     use.dust : should dust (E[B-V]) be fit
##
##
## value
##       coeffs : vector of [distance mod,ebv,peak-to-peak g amp,phase]
ComputeCoeffsPhase <- function(lc,omega,phi,tem,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    m <- m - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    coeffs <- AmpMuDustUpdate(phi,omega,m,t,dust,weights,nb,tem$template_funcs,use.errors,use.dust)
    return(c(coeffs,phi))
}


## make predictions at a set of times (in single band)
##
##
## arguments
##           times : times to make predictions
##            band : band to make prediction
##           omega : frequency
##          coeffs : coefficients from model fit
##             tem : input templates
##
##
## value
##              m  : vector of predicted magnitudes
PredictSingleBand <- function(times,band,omega,coeffs,tem){
    t_temp <- (times*omega + coeffs[4]) %% 1.0
    m <- (coeffs[1] + tem$abs_mag(1/omega,tem)[1,band] + coeffs[2]*tem$dust[band] +
          coeffs[3]*tem$template_funcs[[band]](t_temp))
    return(m)
}

## make predictions at a set of times (in all bands)
##
##
## arguments
##           times : times to make predictions
##           omega : frequency
##          coeffs : coefficients from model fit
##             tem : input templates
##
##
## value
##              m  : matrix of predictions, columns are bands, rows times
PredictAllBand <- function(times,omega,coeffs,tem){
    bands <- colnames(tem$betas)
    m <- vapply(bands,function(x){PredictSingleBand(times,x,omega,coeffs,tem)},rep(0,length(times)))
    return(m)
}


## make predictions at a set of times in particular bands
##
##
## arguments
##           times : times to make predictions
##           bands : bands[ii] is band of observation times[ii]
##           omega : frequency
##          coeffs : coefficients from model fit
##             tem : input templates
##
##
## value
##              m  : vector of predicted magnitudes
PredictTimeBand <- function(times,bands,omega,coeffs,tem){
    m <- rep(0,length(times))
    bands_unique <- unique(bands)
    for(ii in bands_unique){
        m[bands==ii] <- PredictSingleBand(times[bands==ii],ii,omega,coeffs,tem)
    }
    return(m)
}

## grid search across phase at fixed omega.
## possible uses:
## 1. optimize parameter fits AFTER selecting frequency
## 2. test that newton algorithm is finding best parameter fits
##
## arguments
##           lc : light curve, data frame with columns time, band, mag, error
##        omega : frequency
##          tem : input templates
##         phis : grid of phases to try
##   use.errors : should photometric errors be used
##     use.dust : should dust (E[B-V]) be fit
##
##
## value
##          rss : the residual sum of squares at each phase in grid
ComputeRSSPhase <- function(lc,omega,tem,phis=(1:100)/100,use.errors=TRUE,use.dust=TRUE){
    if(use.dust){
        use.dust <- CheckNumberBands(lc)
    }
    CheckLC(lc)
    tem <- CheckTemLC(tem,lc)
    dat <- AugmentData(lc,tem,use.errors)
    m <- dat[[1]]$mag
    dust <- dat[[1]]$dust
    t <- dat[[1]]$time
    weights <- 1 / dat[[1]]$error^2
    nb <- dat[[2]]
    m <- m - rep.int(tem$abs_mag(1/omega,tem)[1,],nb)
    rss_max <- sum((lm(m~dust,weights=weights)$residuals^2)*weights)
    rss <- rep(0,length(phis))
    for(ii in 1:length(phis)){
        coeffs <- AmpMuDustUpdate(phis[ii],omega,m,t,dust,weights,nb,
                                  tem$template_funcs,use.errors,use.dust)
        gammaf <- ConstructGamma(t,nb,phis[ii],omega,tem$template_funcs)
        resid <- m - coeffs[1] - coeffs[2]*dust - coeffs[3]*gammaf
        rss[ii] <- min(sum(weights*resid^2),rss_max)
    }
    return(rss)
}



##### NOTE:
##### below are mostly helper functions that
##### are unlikely to be useful for direct calling

## coerces lc into form for model to fit
AugmentData <- function(lc,tem,use.errors){
    lc <- lc[order(lc$band),]
    nb <- table(lc$band)
    lc$dust <- rep.int(tem$dust,nb)
    ##lc$mag <- lc$mag - rep.int(tem$betas,nb)
    lc$band <- NULL
    if(use.errors){
        lc$error <- sqrt(lc$error^2 + rep.int(tem$model_error,nb)^2) ## adds model error to photometric error
    } else {
        lc$error <- 1
    }
    return(list(lc=lc,nb=nb))
}

## makes \gamma_b(omega*t + \phi) for times t where
## t is ordered (by band) set of times
ConstructGamma <- function(t,nb,phi,omega,temp_funcs){
    t <- (t*omega + phi) %% 1
    bix <- c(0,cumsum(nb))
    gammaf <- rep(0,length(t))
    for(jj in 1:(length(bix)-1)){
        ix1 <- (bix[jj]+1)
        ix2 <- bix[jj+1]
        if(ix2 >= ix1){
            gammaf[ix1:ix2] <-  temp_funcs[[jj]](t[ix1:ix2])
        }
    }
    return(gammaf)
}

## computes a newton update for the (mu,a,d,phi) parameter vector
NewtonUpdate <- function(phi,omega,m,t,dust,weights,nb,template_funcs,templated_funcs,use.errors,use.dust){
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    if(use.dust){
        est <- ComputeBeta(m,dust,gammaf,weights,use.errors)
        mu <- est["mu"]
        a <- est["a"]
        d <- est["d"]
    } else {
        est <- ComputeBetaOne(m,gammaf,weights,use.errors)
        mu <- est["mu"]
        a <- est["a"]
        d <- 0
    }        
    if(a > 0){
        gammafd <- ConstructGamma(t,nb,phi,omega,templated_funcs)
        mp <- m - mu - d*dust
        ## TODO: describe this optimization method somewhere
        ## see FDA book on modified Newton method for phase registration
        if(use.errors){
            del <- sum(gammafd*(mp-a*gammaf)*weights)
            h <- a*sum(gammafd*gammafd*weights)
        } else {
            del <- sum(gammafd*(mp-a*gammaf))
            h <- a*sum(gammafd*gammafd)
        }            
        phi <- (phi + h^{-1}*del) %% 1
    } else {
        a <- 0
        phi <- runif(1)
    }
    out <- c(mu,d,a,phi=phi)
    names(out) <- NULL
    return(out)
}

## update for the (mu,a,d) parameter vector (closed form because phi fixed)
## TODO: can NewtonUpdate just call this function?
AmpMuDustUpdate <- function(phi,omega,m,t,dust,weights,nb,template_funcs,use.errors,use.dust){
    gammaf <- ConstructGamma(t,nb,phi,omega,template_funcs)
    if(use.dust){
        est <- ComputeBeta(m,dust,gammaf,weights,use.errors)
        mu <- est["mu"]
        a <- est["a"]
        d <- est["d"]
    } else {
        est <- ComputeBetaOne(m,gammaf,weights,use.errors)
        mu <- est["mu"]
        a <- est["a"]
        d <- 0
    }        
    if(a < 0) {
        a <- 0
    }
    out <- c(mu,d,a)
    names(out) <- NULL
    return(out)
}

## finds best fitting mu (distance mod) , d (i.e. ebv), a (amplitude)
ComputeBeta <- function(m,dust,gammaf,weights,use.errors){
    X <- cbind(mu=1,d=dust,a=gammaf)
    if (use.errors) {
        B <- t(X)%*%(X*weights)
        d <- t(X)%*%(m*weights)
    } else {
        B <- t(X)%*%X
        d <- t(X)%*%m
    }
    z <- solve(B,d)
    return(z[,1])
}


## finds best fitting mu (distance mod) ,  a (amplitude)
## does not fit dust, so ebv set to 0
ComputeBetaOne <- function(m,gammaf,weights,use.errors){
    X <- cbind(mu=1,a=gammaf)
    if (use.errors) {
        B <- t(X)%*%(X*weights)
        d <- t(X)%*%(m*weights)
    } else {
        B <- t(X)%*%X
        d <- t(X)%*%m
    }
    z <- solve(B,d)
    return(z[,1])
}

## check and make tem and lc consistent
## if lc has bands not in tem, stop
## if lc has fewer bands then tem, get rid
##    of these bands in tem
CheckTemLC <- function(tem,lc){
    if(prod(unique(lc$band) %in% names(tem$dust)) != 1){
        print("template bands are:")
        print(names(tem$dust))
        print("lc is:")
        print(lc)
        stop("all lc bands must match template names")
    }
    bs <- names(tem$dust)[names(tem$dust) %in% unique(lc$band)]
    if(length(bs) < length(tem$dust)){
        tem$betas <- tem$betas[,bs,drop=FALSE]
        tem$dust <- tem$dust[bs]
        tem$model_error <- tem$model_error[bs]
        tem$templates <- tem$templates[bs,]
        tem$templatesd <- tem$templatesd[bs,]
        tem$template_funcs <- tem$template_funcs[bs]
        tem$templated_funcs <- tem$templated_funcs[bs]
    }
    return(tem)
}

## useful for calling these functions from python, see template.py for example
TBMEtoLC <- function(time,band,mag,error){
    return(data.frame(time,band,mag,error,stringsAsFactors=FALSE))
}

## if lc has only one band, need to use ComputeBetaOne rather than ComputeBeta
## this function checks / warns if using single band
CheckNumberBands <- function(lc){
    if(length(unique(lc$band))==1){
        print("warning: light curve has only 1 band, setting E[B-V] = 0 (i.e. assume no dust). the distance modulus is now the band mean and has no physical interpretation unless the light curve was already dust corrected. set use.dust=FALSE to prevent this warning message from being displayed.")
        use.dust <- FALSE
    } else {
        use.dust <- TRUE
    }
    return(use.dust)
}

## checks that lc has correct column structure / sufficient rows
CheckLC <- function(lc){
    if(nrow(lc)<5){
        stop("lc must have at least 5 rows (observations)")
    }
    if(sum(names(lc) %in% c("time","band","mag","error")) != 4){
        stop("lc must have columns named: time, band, mag, error")
    }
    if(typeof(lc$time)!="double"){
        stop("typeof(lc$time) must be double")
    }
    if(typeof(lc$mag)!="double"){
        stop("typeof(lc$mag) must be double")
    }
    if(typeof(lc$error)!="double"){
        stop("typeof(lc$error) must be double")
    }
    if(typeof(lc$band)!="character"){
        stop("typeof(lc$band) must be character")
    }
}

    
