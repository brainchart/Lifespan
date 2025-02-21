source("100.common-variables.r")
source("101.common-functions.r")

source("300.variables.r")
source("301.functions.r")

FIT <- readRDS("Share/FIT_GMV.rds")

POP.CURVE.LIST <- list(AgeTransformed=seq(log(90),log(365*95),length.out=2^4),sex=c("Female","Male"))
POP.CURVE.RAW <- do.call( what=expand.grid, args=POP.CURVE.LIST )

CURVE <- Apply.Param(NEWData=POP.CURVE.RAW, FITParam=FIT$param )

plot( PRED.m500.pop ~ AgeTransformed, data=CURVE[CURVE$sex=="Female",],type="l")

plot( I(10000*PRED.m500.pop) ~ AgeTransformed, data=CURVE[CURVE$sex=="Female",],type="l")


