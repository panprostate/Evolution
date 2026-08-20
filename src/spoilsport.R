# Helper function from "https://github.com/morrislab/SpoilSport" to perform power correction of number of mutations

library(bbmle)
library(emdbook)


wcc <- function(purity,ploidy,rd,T,props){
    T = T-1
    out = data.frame(
                    sc_n=integer(),
                    cp=double(),
                    rd = integer(),
                    ploidy = double(),
                    sf = double(),
                    tp = double(), stringsAsFactors=FALSE)

    for (i in 1:length(props)){
        bp = props[i]/(purity*ploidy+(1.0-purity)*2.0)
        trunc_fit <- function(tbp){
            mu=0
            spdf = 0
            for (j in (T+1):rd){ 
                p = dbinom(j,rd,tbp)
                mu = mu + p*j
                spdf= spdf + p
            }  
            abs(bp - (mu/spdf)/rd) + (1-tbp)*.4
            }
        m0 = optim(0.3,trunc_fit,method="Brent",upper=c(1.0),lower=c(0.01))
        tbp = m0$par

        sp = 0
        for (j in 0:T){
            sp = sp + dbinom(j,rd,tbp)
        }
        sf = 1/(1-sp)
        tbp = tbp*(purity*ploidy+(1.0-purity)*2.0)
        if (i == 1){
            sf = 1
            tbp = props[i]
        }
        out = rbind(out,data.frame(i,props[i],sf,tbp))

    }
    colnames(out) = c("sc_n","cp","sf","tp")
    out
}