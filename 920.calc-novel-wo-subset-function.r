# rm(list=ls())

if( 0 ) {
    dir("RDS")
    NEWDATA <- read.csv("RDS/dummydata.csv")
    NEWDATA <- read.csv("RDS/saashi_abide_dummdata_test_WMV.csv")
    FIT <- readRDS("RDS/Model/FIT_WMV.rds")
    BOOT <- readRDS("RDS/Model/BOOT_WMV.rds")

    RESULT <- Calc.Novel( NewData=NEWDATA, Fit=FIT, BootFit=BOOT[1:20], Apply=TRUE, NumberCores=6 )
    NewData=NEWDATA; Fit=FIT; BootFit=BOOT[1:20]; Apply=TRUE; NumberCores=6
}

if( 0 ) {
    NewData.all <- read.csv( file="TEST/NEW/NEWDATA.csv" )
    NewData.all$meanCT2 <- NewData.all$WMV/100 ## for testing only

    NewData.v1 <- NewData.all[1:10,] ## for testing
    NewData.v1$fs_version <- "FS-TEST"
    NewData.v1$fs_version[1:3] <- NA ## for testing
    NewData.v1$sex[1:3] <- NA ## for testing

    Fit.GMV <- readRDS( file="TEST/GMV/FIT.EXTRACT.rds" )
    Boot.GMV <- readRDS( file="TEST/GMV/BOOT.EXTRACT.rds" )

    Fit.CT <- readRDS( file="TEST/CT/FIT.EXTRACT.rds" )


    RESULT <- Calc.Novel( NewData.v1, Fit.CT, Apply=TRUE, NumberCores=6 )

    RESULT <- Calc.Novel( NewData.all, Fit=Fit.GMV, BootFit=Boot.GMV, Apply=TRUE, NumberCores=6 )

    RESULTA <- Calc.Novel( NewData.all, Fit=Fit.GMV, BootFit=NULL, Apply=TRUE )
    RESULT1 <- Calc.Novel( NewData.all[NewData.all$study=="DUMMY1",], Fit=Fit.GMV, BootFit=Boot.GMV, Apply=TRUE )
    RESULT2 <- Calc.Novel( NewData.all[NewData.all$study=="DUMMY2",], Fit=Fit.GMV, BootFit=NULL, Apply=TRUE )

    RESULTA$fit$param$mu$ranef[c("DUMMY1","DUMMY2")]
    RESULT1$fit$param$mu$ranef[c("DUMMY1","DUMMY2")]
    RESULT2$fit$param$mu$ranef[c("DUMMY1","DUMMY2")]
    head(RESULTA$data)
}

## NOTE: remove Cond.Base argument, use "level" attribute

Calc.Novel <- function( NewData, Fit, BootFit=NULL, Apply=TRUE, NumberCores=1 ) {

    require("gamlss") ## need gamlss.family objects in all scripts

    if( NumberCores>1 ) {
        require("parallel") ## start with the simply mclapply (multicore-lapply from the parallel package)
    }
    
    ## ===== SECTION COPIED FROM 102 GAMLSS RECODE =====
    ##
    ##
    bfpNA <- function (x, powers = c(1, 2), shift = 0, scale = 1) {
        nobs <- length(x)
        npoly <- length(powers)
        X <- matrix(0, nrow = nobs, ncol = npoly)
        if (is.null(scale) | is.null(shift)) {
            stop("WARNING: Using automatic scale/shift will invalidate future refitting")
            out <- fp.scale(x)
            shift <- out$shift
            scale <- out$scale
        }
        x1 <- ifelse(powers[1] != rep(0, nobs), x^powers[1], log(x))
        X[, 1] <- x1
        if (npoly >= 2) {
            for (i in 2:npoly) {
                if (powers[i] == powers[(i - 1)]) 
                    x2 <- log(x) * x1
                else x2 <- ifelse(powers[i] != rep(0, nobs), x^powers[i], 
                                  log(x))
                X[, i] <- x2
                x1 <- x2
            }
        }
        X
    }
    bfp <- function( ... ) {stop("Default bfp() function cannot handle NAs. We have masked with this fatal error. Use bfpNA() instead. ")}
    dGGalt <- dGG
    pGGalt <- pGG
    qGGalt <- qGG
    rGGalt <- rGG
    GGalt <- function (mu.link = "log", sigma.link = "log", nu.link = "identity") {
        mstats <- checklink("mu.link", "GGalt", substitute(mu.link), c("1/mu^2", "log", "identity"))
        dstats <- checklink("sigma.link", "GGalt", substitute(sigma.link), c("inverse", "log", "identity"))
        vstats <- checklink("nu.link", "GGalt", substitute(nu.link), c("1/nu^2", "log", "identity"))
        GGalt.mean <- function (mu, sigma, nu) {
            TOP <- log(mu) + lgamma(1/(sigma^2 * nu^2) + 1/nu)
            BOTTOM <-  (-1)*((1/nu)*log((sigma^2 * nu^2))) + lgamma(1/(sigma^2 * nu^2))
            ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 1), exp(TOP-BOTTOM), Inf)
        }
        GGalt.variance <- function (mu, sigma, nu) {
            AA <- 2*log(mu) - ((2/nu)*(-1)*log( (sigma^2 * nu^2) )) - 2*lgamma(1/(sigma^2 * nu^2))
            ww <- lgamma(1/(sigma^2 * nu^2) + 2/nu) + lgamma(1/(sigma^2 * nu^2))
            uu <- 2*lgamma(1/(sigma^2 * nu^2) + 1/nu)
            BB <- ww + log( (1 - exp( uu - ww )) )
            YES <- AA + BB
            ifelse(nu > 0 | (nu < 0 & sigma^2 * abs(nu) < 0.5), ifelse(is.nan(BB),NA,exp( YES )), Inf)
        }

        structure(list(family = c("GGalt", "generalised Gamma Lopatatsidis-Green (altered)"),
                       parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                       nopar = 3,
                       type = "Continuous",
                       ##
                       mu.link = as.character(substitute(mu.link)), 
                       sigma.link = as.character(substitute(sigma.link)),
                       nu.link = as.character(substitute(nu.link)), 
                       ##
                       mu.linkfun = mstats$linkfun,
                       sigma.linkfun = dstats$linkfun, 
                       nu.linkfun = vstats$linkfun,
                       ##
                       mu.linkinv = mstats$linkinv, 
                       sigma.linkinv = dstats$linkinv,
                       nu.linkinv = vstats$linkinv,
                       ##
                       mu.dr = mstats$mu.eta,
                       sigma.dr = dstats$mu.eta,
                       nu.dr = vstats$mu.eta,
                       ##
                       dldm = function(y, mu, sigma, nu) {
                           z <- (y/mu)^nu
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           dldm <- ifelse(abs(nu) > 1e-06, (z - 1) * theta * nu/mu, (1/(mu * (sigma^2)) * (log(y) - log(mu))))
                           dldm },
                       d2ldm2 = function(mu, sigma, nu) {
                           d2ldm2 <- ifelse(abs(nu) > 1e-06, -1/((mu^2) * (sigma^2)), -(1/(mu^2 * sigma^2)))
                           d2ldm2 },
                       dldd = function(y, mu, sigma, nu) {
                           z <- (y/mu)^nu
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           dldd <- ifelse(abs(nu) > 1e-06, -2 * theta * (log(theta) + 1 + log(z) - z - digamma(theta))/sigma, -(1/sigma) + (1/sigma^3) * (log(y) - log(mu))^2)
                           dldd },
                       d2ldd2 = function(y, mu, sigma, nu) {
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           d2ldd2 <- ifelse(abs(nu) > 1e-06, 4 * (theta/(sigma^2)) * (1 - theta * trigamma(theta)), -2/sigma^2)
                           d2ldd2 },
                       dldv = function(y, mu, sigma, nu) {
                           z <- (y/mu)^nu
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           dldv <- (1/nu) * (1 + 2 * theta * (digamma(theta) + z - log(theta) - 1 - ((z + 1)/2) * log(z)))
                           dldv },
                       d2ldv2 = function(y, mu, sigma, nu) {
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           d2ldv2 <- -(theta/nu^2) * (trigamma(theta) * (1 + 4 * theta) - (4 + 3/theta) - log(theta) * (2/theta - log(theta)) + digamma(theta) * (digamma(theta) + (2/theta) - 2 * log(theta)))
                           d2ldv2 },
                       d2ldmdd = function(y) rep(0, length(y)),
                       d2ldmdv = function(y, mu, sigma, nu) {
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           ddd <- (theta/mu) * (digamma(theta) + (1/theta) - log(theta))
                           ddd },
                       d2ldddv = function(y, mu, sigma, nu) {
                           theta <- 1/(sigma^2 * abs(nu)^2)
                           d2ldddv <- -2 * sign(nu) * theta^(3/2) * (2 * theta * trigamma(theta) - (1/theta) - 2)
                           d2ldddv },
                       G.dev.incr = function(y, mu, sigma, nu, ...) -2 * dGGalt(y, mu = mu, sigma = sigma, nu = nu, log = TRUE), ## EDIT
                       ##
                       rqres = expression(rqres(pfun = "pGGalt", type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu)), ## EDIT
                       ##
                       mu.initial = expression(mu <- (y + mean(y))/2),
                       sigma.initial = expression(sigma <- rep(1, length(y))),
                       nu.initial = expression(nu <- rep(1, length(y))), mu.valid = function(mu) all(mu > 0),
                       ##
                       sigma.valid = function(sigma) all(sigma > 0),
                       nu.valid = function(nu) TRUE, 
                       y.valid = function(y) all(y > 0),
                       ##
                       mean = GGalt.mean, ## EDIT
                       variance = GGalt.variance ## EDIT
                       ##
                       ),
                  class = c("gamlss.family", "family")
                  )
    }

    
    if( missing(NewData) ) stop("Must supply NewData argument")
    if( missing(Fit) ) stop("Must supply a Fit object")
    if( ! any(class(Fit$param)=="ParamObj") ) stop("Fit object does not contain a ParamObj element")
    if( !is.null(BootFit) ) {
        WHICH.NULL <- which(sapply(BootFit,function(Z){is.null(Z$param)}))
        if( length( WHICH.NULL ) > 0 ) {
            warning("Dropping elements of BootFit without a converged ParamObj: ",length(WHICH.NULL)," elements dropped (",paste0(WHICH.NULL,collapse=","),")")
            BootFit <- BootFit[ -WHICH.NULL ]
        } 
        if( !all(sapply(BootFit,function(X){any(class(X$param)=="ParamObj")})) ) stop("BootFit object does not contain a list of ParamObj elements")
    }

    if( all(is.na(Fit$param$mu$ranef)) & all(is.na(Fit$param$sigma$ranef)) ) stop("No random-effects in mu-component or sigma-component")
    

    lCOND <- attr(Fit$param,"model")$covariates$COND
    lBY <- attr(Fit$param,"model")$covariates$BY
    lOTHER <- attr(Fit$param,"model")$covariates$OTHER
    lRANEF <- attr(Fit$param,"model")$covariates$RANEF
    if( length(lRANEF) != 1 ) stop("This function assumes a single random-effect covariate (that may appear in each component)")
    if( ! all(c(lBY,lCOND,lRANEF,lOTHER) %in% names(NewData)) ) stop("NewData missing covariate and/or random-effect columns")

    for( LAB in names(attr(Fit$param,"levels")) ) {
        NewData[,sprintf("%s.original",LAB)] <- NewData[,LAB]
        NewData[,LAB] <- factor( as.character(NewData[,LAB]), levels=attr(Fit$param,"levels")[[LAB]] )
    }
    
    for( LAB in names(attr(Fit$param,"transformations")) ) {
        if( !(attr(Fit$param,"transformations")[[LAB]]$OriginalName %in% names(NewData) ) ) {stop("Source column for transformation not in NewData")}
        Local.Function <- attr(Fit$param,"transformations")[[LAB]]$toTransformed
        NewData[,attr(Fit$param,"transformations")[[LAB]]$TransformedName] <- Local.Function( NewData[,attr(Fit$param,"transformations")[[LAB]]$OriginalName] )

    }
    
    
    lY <- attr(Fit$param,"model")$covariates$Y
    lX <- attr(Fit$param,"model")$covariates$X
    if( ! all(c(lY) %in% names(NewData)) ) stop("NewData missing outcome column")
    if( ! all(c(lX) %in% names(NewData)) ) stop("NewData missing covariate and/or random-effect columns")
    
    Cond.Base <- attr(Fit$param,"levels")[[lCOND]]
    if( (length(Cond.Base)!=1) || (!is.character(Cond.Base)) ) stop("Ambiguous condition baseline detected")
    if( ! Cond.Base %in% levels(NewData[,lCOND]) ) stop("NewData does not include any indiviudals with ",lCOND," equal to ",Cond.Base)


    for( LAB in names(attr(Fit$param,"levels")) ) {
        if( LAB == lCOND ) { next } ## must skip the condition column, since witihn LOCAL there is only one level (healthy controls)
        if( !is.factor( NewData[,LAB] ) ) {stop("NewData column is not a factor when Fit object expects it to be")}
        
        if( LAB %in% names(attr(Fit$param,"contrasts")) ) {
            contrasts(NewData[,LAB]) <- attr(Fit$param,"contrasts")[[LAB]]
        } else if( is.ordered( NewData[,LAB] ) ) {
            contrasts(NewData[,LAB]) <- "contr.poly"
        } else {
            contrasts(NewData[,LAB]) <- "contr.treatment"
        }
    }
    attr(NewData,"na.action") <- "na.pass"

    
    LOCAL.ALL <- NewData[ which(NewData[,lCOND]==Cond.Base), c(lY,lX,lBY,lCOND,lRANEF,lOTHER)]

    if( NROW(LOCAL.ALL)==0 ) {stop("There are no control observations on which to derive a novel estimate")}
    
    ##
    ## Code assumes nested random-effect structure, mu contains sigma contains nu contains tau
    ## Hence, we need to only consider mu-component random-effects for missing levels
    if( is.factor(LOCAL.ALL[,lRANEF]) ) {
        NewLevels <- levels(LOCAL.ALL[,lRANEF])[! levels(LOCAL.ALL[,lRANEF]) %in% names(Fit$param$mu$ranef)]
    } else if( is.character(LOCAL.ALL[,lRANEF]) ) {
        NewLevels <- unique(LOCAL.ALL[,lRANEF])[ ! unique(LOCAL.ALL[,lRANEF]) %in% names(Fit$param$mu$ranef)]
    } else {
        stop("Code assumes RANEF is a factor or character column")
    }

    LOCAL <- LOCAL.ALL[ which(LOCAL.ALL[,lRANEF] %in% NewLevels), ]


    LOCAL.MISSING <- which( (rowSums(is.na(LOCAL[ , c(lY,lX,lBY,lRANEF)]))>0) ) ## lCOND and lOTHER are excluded
    if( length(LOCAL.MISSING)>0 ) {
        warning("There are missing values in the key columns: ",paste0(c(lY,lX,lBY,lRANEF),collapse=", "),". Dropping ",length(LOCAL.MISSING)," rows from ML estimation.")
        LOCAL <- LOCAL[-LOCAL.MISSING,]
    }
    attr(LOCAL,"na.action") <- "na.pass"    
    
    Make.Known <- function( lFIT, lLOCAL, lOUT, lRANEF ) {
        ## lFIT=Fit; lLOCAL=LOCAL; lOUT=lY
        
        KNOWN <- list()
        KNOWN$Y <- lLOCAL[,lOUT]
        attr(KNOWN,"NewLevels") <- NewLevels
        attr(KNOWN,"family") <- lFIT$param$family
        FamilyComponents <- names(get(lFIT$param$family)()$parameters)
        for( IDX in 1:length(FamilyComponents) ) {
            LAB <- FamilyComponents[IDX]
            MM <- model.matrix(as.formula(paste("~",lFIT$param[[LAB]]$equ$fixef)), data=lLOCAL )
            for( REPLACE.NAME in names(attr(lFIT$param,"model")$contrasts) ) {
                if( attr(lFIT$param,"model")$contrasts[[REPLACE.NAME]] == "contr.sum" ) {
                    ## if all NA, replace with zero
                    COMPONENTS <- grep(REPLACE.NAME,colnames(MM),value=TRUE)
                    WHICH <- apply(MM[,COMPONENTS],1,function(x){all(is.na(x))})
                    MM[WHICH,COMPONENTS] <- 0
                } else {
                    stop("NewData is missing entries for ",REPLACE.NAME," (factor column) and we do not know what to set them as?")
                }
            }
            
            BETA <- matrix(lFIT$param[[LAB]]$fixef,ncol=1,dimnames=list(names(lFIT$param[[LAB]]$fixef),NULL))
            KNOWN[[LAB]]  <- list(fixef=MM %*% BETA[match(colnames(MM),rownames(BETA))],
                                  sigma=if( is.na(lFIT$param[[LAB]]$equ$ranef) ) {NA} else { lFIT$param[[LAB]]$sigma },
                                  map=match(lLOCAL[,lRANEF],NewLevels) )
            if( !is.na(lFIT$param[[LAB]]$equ$ranef) ) {
                KNOWN[[LAB]]$theta <-  1:length(unique(KNOWN[[LAB]]$map)) + if(IDX>1){ sum(sapply(KNOWN[FamilyComponents[1:(IDX-1)]],function(Z){length(Z$theta)})) }else{0}
            } else {
                KNOWN[[LAB]]$theta <- NULL
            }
        }
        attr(KNOWN,"FamilyComponents") <- FamilyComponents
        attr(KNOWN,"ThetaComponents") <- FamilyComponents[sapply(KNOWN[FamilyComponents],function(X){!is.null(X$theta)})]
        attr(KNOWN,"FixefComponents") <- setdiff(attr(KNOWN,"FamilyComponents"),attr(KNOWN,"ThetaComponents"))
        return(KNOWN)
    }    

    Make.Theta <- function (lknown ) {
        THETA <- rep( 0, length.out=sum(lengths(lapply( lknown[attr(lknown,"FamilyComponents")], function(X){X$theta} ))) )
        return(THETA)
    }

    MLE.Func <- function( theta, Known, Return=c("list","value","optim")[3] ) {
        LL <- list()
        
        for( LAB in attr(Known,"ThetaComponents") ) {
            LL[[paste0(LAB,".sigma")]] <- dnorm(x=theta[ Known[[LAB]]$theta ],mean=0,sd=Known[[LAB]]$sigma,log=TRUE) ## this can be a vector of length>1
        }

        NOVEL <- list()
        for( LAB in attr(Known,"ThetaComponents") ) {
            NOVEL[[LAB]] <- Known[[LAB]]$fixef + theta[ Known[[LAB]]$theta ][ Known[[LAB]]$map ]
        }
        for( LAB in attr(Known,"FixefComponents") ) {
            NOVEL[[LAB]] <- Known[[LAB]]$fixef
        }
        
        DIST <- get( attr(Known,"family") )()
        lARGS <- list()
        CHECK <- DIST$parameters
        for( LAB in attr(Known,"FamilyComponents") ) {
            lARGS[[LAB]] <- DIST[[sprintf("%s.linkinv",LAB)]]( NOVEL[[LAB]] )
            CHECK[[LAB]] <- DIST[[sprintf("%s.valid",LAB)]]( lARGS[[LAB]] )
        }
        if( !all(unlist(CHECK)) ) {stop("Failed distribution parameter checks")}

        LL$out <- do.call( what=get(paste0("d",attr(Known,"family"))), args=c(lARGS,list(x=Known$Y),log=TRUE))
        if(Return=="list"){
            return(LL)
        } else if (Return=="optim") {
            return( -1*sum(unlist(LL)) )            
        } else if (Return=="value") {
            return( sum(unlist(LL)) )
        } else {
            stop("Invalid return mode")
        }
    }

    Optim.Novel <- function( lFIT, lLOCAL, lOUT ) {
        ## lFIT=Fit; lLOCAL=LOCAL; lOUT=lY
        
        ## This function assumes that the following objects are in scope: LOCAL, lY
        ## (should probably re-write these to explicitly pass these arguments, but since this is wrapping all calc-novel steps into a single function)
        localKnown <- Make.Known( lFIT=lFIT, lLOCAL=lLOCAL, lOUT=lOUT, lRANEF=lRANEF )
        localTheta <- Make.Theta( localKnown )
        localOpt <- optim(par=localTheta, fn=MLE.Func, Known=localKnown,
                          method=if(length(localTheta)==1){"Brent"}else{"Nelder-Mead"},
                          lower=if(length(localTheta)==1){-1000}else{-Inf},
                          upper=if(length(localTheta)==1){ 1000}else{ Inf})

        if( (length(localTheta)>1) & (localOpt$convergence > 0) ) {
            ## need to try alternative method (for longer length Theta vectors
            ## (Also increase number of iterations)
            cat('WARNING: no MLE convergence using Brent or Nelder-Mead, reverting to Broyden–Fletcher–Goldfarb–Shanno optimisation \n')
            localOpt  <- optim(par=localTheta, fn=MLE.Func, Known=localKnown, control=list(maxit=1000), method="BFGS" )
        }
        
        if( localOpt$convergence > 0 ) {
            cat("Optim.Novel() did not converge for (",if(!is.null(lFIT$offset)){lFIT$offset}else{"main"},")","\n")
            stop("Currently the code that applies the bootstrap parameters to the newdata is not safe given failed optim() convergence. Hence we stop.")
            return(lFIT)
        } else {
            EXPANDED <- lFIT
            for( LAB in attr(localKnown,"ThetaComponents") ) {
                names(EXPANDED$param[[LAB]]$ranef.TYPE) <- names(EXPANDED$param[[LAB]]$ranef) ## fix legacy issue of not naming the ranef.TYPE vector
                EXPANDED$param[[LAB]]$ranef[ attr(localKnown,"NewLevels")[unique(localKnown[[LAB]]$map)] ] <- localOpt$par[ localKnown[[LAB]]$theta ]
                EXPANDED$param[[LAB]]$ranef.TYPE[ attr(localKnown,"NewLevels")[unique(localKnown[[LAB]]$map)] ] <- "novel"
            }
            return(EXPANDED)
        }
    }

    cat("NewData inspected and any novel studies identified. Finding Maximum Likelihood estimates for random-effects. \n")
    
    if( length(NewLevels)> 0 ) {

        if( !is.null(BootFit) ) {
            if( NumberCores > 1 ) {
                cat(sprintf("Computing novel random-effects for bootstrap replicates using %i cores. This may take some time.\n",NumberCores))
                TEMP.BOOT <- mclapply( X=BootFit, FUN=Optim.Novel, lLOCAL=LOCAL, lOUT=lY, mc.cores=NumberCores )
                
            } else {
                warning("Estimating novel study random-effects for bootstrap replicates in serial, this may be very slow.\n(Consider using parallelisation option)")

                TEMP.BOOT <- lapply( BootFit, Optim.Novel, lLOCAL=LOCAL, lOUT=lY )
            }
            cat("Novel bootstrap estimaton complete.\n")
        } else {
            TEMP.BOOT <- NULL
        }
        
        RETURN <- list(fit=Optim.Novel( lFIT=Fit, lLOCAL=LOCAL, lOUT=lY ),
                       boot=TEMP.BOOT
                       )

    } else {
        ## return input unchanged if no unknown levels
        RETURN <- list(fit=( Fit ),
                       boot=( BootFit )
                       )

        
    }

    cat("Novel study random-effects calculated.\n")
    
    if( Apply==TRUE ) {
        
        NewData[,sprintf("%s.q.wre",lY)] <- NA_real_
        NewData[,sprintf("%s.normalised",lY)] <- NA_real_
        NewData[,sprintf("%s.m500.wre",lY)] <- NA_real_
        NewData[,sprintf("%s.m500.pop",lY)] <- NA_real_

        WHICH <- which( (rowSums(is.na(NewData[ , c(lY,lX,lBY,lRANEF)]))==0) ) ## lCOND and lOTHER are excluded
        if( length(WHICH) != NROW(NewData) ) {warning("There are missing values in the key columns: ",paste0(c(lY,lX,lBY,lRANEF),collapse=", "),". Skipping ",NROW(NewData)-length(WHICH)," rows for ")}


        ApplyParam <- function( X, APPLYDF,ALL=FALSE ) {
            FamilyComponents <- names(get(X$param$family)()$parameters)
            ARGS <- list()
            POPS <- list()
            for( LAB in FamilyComponents ) {
                Model.Formula <- as.formula(paste("~",X$param[[LAB]]$equ$fixef))
                Model.Matrix <- model.matrix( object=Model.Formula, data=APPLYDF )
                for( REPLACE.NAME in names(attr(X$param,"model")$contrasts) ) {
                    if( attr(X$param,"model")$contrasts[[REPLACE.NAME]] == "contr.sum" ) {
                        ## if all NA, replace with zero
                        COMPONENTS <- grep(REPLACE.NAME,colnames(Model.Matrix),value=TRUE)
                        REPLACE.WHICH <- apply(Model.Matrix[,COMPONENTS],1,function(x){all(is.na(x))})
                        Model.Matrix[REPLACE.WHICH,COMPONENTS] <- 0
                    } else {
                        stop("NewData is missing entries for ",REPLACE.NAME," (factor column) and we do not know what to set them as?")
                    }
                }
            
                BETA <- matrix(X$param[[LAB]]$fixef,ncol=1,dimnames=list(names(X$param[[LAB]]$fixef),NULL))
                APPLYDF[,sprintf("%s.pop",LAB)] <- as.numeric(Model.Matrix %*% BETA[match(colnames(Model.Matrix),rownames(BETA))])

                if( !is.na( X$param[[LAB]]$equ$ranef ) ) {
                    APPLYDF[,sprintf("%s.ranef",LAB)] <- X$param[[LAB]]$ranef[ as.character(APPLYDF[,lRANEF]) ]
                    APPLYDF[,sprintf("%s.wre",LAB)] <- APPLYDF[,sprintf("%s.pop",LAB)] + APPLYDF[,sprintf("%s.ranef",LAB)]
                } else {
                    APPLYDF[,sprintf("%s.ranef",LAB)] <- Inf
                    APPLYDF[,sprintf("%s.wre",LAB)] <- APPLYDF[,sprintf("%s.pop",LAB)]
                }
                ARGS[[LAB]] <- get(X$param$family)()[[sprintf("%s.linkinv",LAB)]]( APPLYDF[,sprintf("%s.wre",LAB)] )
                POPS[[LAB]] <- get(X$param$family)()[[sprintf("%s.linkinv",LAB)]]( APPLYDF[,sprintf("%s.pop",LAB)] )
            }

            APPLYDF[,sprintf("%s.m500.wre",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(ARGS,list(p=0.500)))
            APPLYDF[,sprintf("%s.m500.pop",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(POPS,list(p=0.500)))
            APPLYDF[,sprintf("%s.l025.wre",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(ARGS,list(p=0.025)))
            APPLYDF[,sprintf("%s.l025.pop",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(POPS,list(p=0.025)))
            APPLYDF[,sprintf("%s.u975.wre",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(ARGS,list(p=0.975)))
            APPLYDF[,sprintf("%s.u975.pop",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(POPS,list(p=0.975)))                        
            APPLYDF[,sprintf("%s.q.wre",lY)] <- do.call(what=get(paste0("p",X$param$family)), args=c(ARGS,list(q=APPLYDF[,lY])))
            APPLYDF[,sprintf("%s.normalised",lY)] <- do.call(what=get(paste0("q",X$param$family)), args=c(POPS,list(p=APPLYDF[,sprintf("%s.q.wre",lY)])))
            
            RCOLUMNS <- c(paste0(lY,c(".m500.wre",".m500.pop",".u975.wre",".u975.pop",".l025.wre",".l025.pop",".q.wre",".normalised")),
                          paste0(FamilyComponents,c(".pop")),paste0(FamilyComponents,c(".wre")))
            if( !ALL ) {
                return( APPLYDF[,RCOLUMNS] )
            } else {
                return(APPLYDF)
            }
        }


        if( NumberCores > 1 ) {
            if(!is.null(BootFit)) cat(sprintf("Computing bootstrap intervals using %i cores. This may take some time.\n",NumberCores))
            TEMP.LIST <- mclapply( c(list(RETURN$fit),RETURN$boot), ApplyParam, APPLYDF=NewData[WHICH,], mc.cores=NumberCores )
        } else {
            if(!is.null(BootFit)) warning("Computing bootstrap intervals in serial, this may be very slow.\n(Consider using parallelisation option)")
            TEMP.LIST <- lapply( c(list(RETURN$fit),RETURN$boot), ApplyParam, APPLYDF=NewData[WHICH,] )
        }
        if(!is.null(BootFit)) cat("Bootstrap intervals complete\n")
        
        RETURN$data <- NewData
        RETURN$data[WHICH,names(TEMP.LIST[[1]])] <- TEMP.LIST[[1]]

        if(!is.null(BootFit)) {
            for( JLAB in paste0(lY,c(".m500.wre",".m500.pop",".u975.wre",".u975.pop",".l025.wre",".l025.pop",".q.wre",".normalised")) ) {
                RETURN$data[,paste0(JLAB,c(".bootl",".bootu"))] <- t(apply( sapply( TEMP.LIST[-1], function(Z) { Z[,JLAB]} ),1,FUN=quantile,probs=c(0.025,0.975)))
            }
        }
        
    }
    
    return(RETURN)
}
