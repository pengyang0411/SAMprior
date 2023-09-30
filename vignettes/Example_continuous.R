## ---- include=FALSE-----------------------------------------------------------
library(SAMprior)
library(knitr)
knitr::opts_chunk$set(
    fig.width = 1.62*4,
    fig.height = 4
    )
## setup up fast sampling when run on CRAN
is_CRAN <- !identical(Sys.getenv("NOT_CRAN"), "true")
## NOTE: for running this vignette locally, please uncomment the
## following line:
## is_CRAN <- FALSE
.user_mc_options <- list()
if (is_CRAN) {
    .user_mc_options <- options(RBesT.MC.warmup=250, RBesT.MC.iter=500, RBesT.MC.chains=2, RBesT.MC.thin=1, RBesT.MC.control=list(adapt_delta=0.9))
}

## ----results="D_h",echo=FALSE-------------------------------------------------
set.seed(123)
std <- function(x) sd(x)/sqrt(length(x))
df_1 <- rnorm(40, 0, 3); 
df_2 <- rnorm(50, 0, 3); 
df_3 <- rnorm(60, 0, 3); 
dat <- data.frame(study = c(1,2,3),
                  n = c(40, 50, 60),
                  mean = round(c(mean(df_1), mean(df_2), mean(df_3)), 3),
                  se = round(c(std(df_1), std(df_2), std(df_3)), 3))
kable(dat)

## ---- message=FALSE-----------------------------------------------------------
sigma = 3
# load R packages
library(ggplot2)
theme_set(theme_bw()) # sets up plotting theme
set.seed(22)
map_mcmc <- gMAP(cbind(mean, se) ~ 1 | study, 
                 weights=n,data=dat,
                 family=gaussian,
                 beta.prior=cbind(0, sigma),
                 tau.dist="HalfNormal",tau.prior=cbind(0,sigma/2))

map_automix <- automixfit(map_mcmc)
map_automix
plot(map_automix)$mix

## ---- message=FALSE-----------------------------------------------------------
set.seed(234)
sigma        <- 3 ## Standard deviation in the current trial
data.crt     <- rnorm(30, mean = 0.4, sd = sigma)
wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.5 * sigma,
                   data = data.crt)
cat('SAM weight: ', wSAM)

## ---- message=FALSE-----------------------------------------------------------
wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.5 * sigma,
                   method.w = 'PPR',
                   prior.odds = 3/7,
                   data = data.crt)
cat('SAM weight: ', wSAM)

## ---- echo=FALSE, message=FALSE, warning=FALSE--------------------------------
weight_grid <- seq(-3, 3, by = 0.3)
weight_res  <- lapply(weight_grid, function(x){
  res <- c()
  for(s in 1:300){
    data.control <- rnorm(n = 35, mean = x, sd = sigma)
    res <- c(res, SAM_weight(if.prior = map_automix,
                             delta = 0.5 * sigma,
                             data = data.control))
    
  }
  mean(res)
})
df_weight <- data.frame(grid   = weight_grid,
                        weight = unlist(weight_res))
qplot(grid, weight, data = df_weight, geom = "line", main= "SAM Weight") +
  xlab('Sample mean from control trial')+ ylab('Weight') +
  geom_vline(xintercept = summary(map_automix)['mean'], linetype = 2, col = 'blue') 

## ---- message=FALSE-----------------------------------------------------------
unit_prior <- mixnorm(nf.prior = c(1, summary(map_automix)['mean'], sigma))
SAM.prior <- SAM_prior(if.prior = map_automix, 
                       nf.prior = unit_prior,
                       weight = wSAM, sigma = sigma)
SAM.prior

## ---- message=FALSE-----------------------------------------------------------
set.seed(123)
# weak_prior <- mixnorm(c(1, summary(map_automix)[1], 1e4))
TypeI <- get_OC(if.prior = map_automix,    ## MAP prior from historical data
                nf.prior = unit_prior,     ## Weak-informative prior for treatment arm
                delta    = 0.5*sigma,      ## CSD for SAM prior
                n        = 35, n.t = 70,   ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,           ## Number of trials simulated
                if.MAP   = TRUE,           ## Output robust MAP prior for comparison
                weight   = 0.5,            ## Weight for robust MAP prior
                ## Mean for control and treatment arms
                theta    = c(0, 0,    -2, 4),
                theta.t  = c(0, -0.1, -2, 4),
                sigma    = sigma
                  )
kable(TypeI)

## ---- message=FALSE-----------------------------------------------------------
set.seed(123)
Power <- get_OC(if.prior = map_automix,    ## MAP prior based on historical data
                nf.prior = unit_prior,     ## Non-informative prior for treatment arm
                delta    = 0.5*sigma,      ## CSD for SAM prior
                n        = 35, n.t = 70,   ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,           ## Number of trials simulated
                if.MAP   = TRUE,           ## Output robust MAP prior for comparison
                weight   = 0.5,            ## Weight for robust MAP prior
                ## Mean for control and treatment arms
                theta    = c(0, 0.1, 0.5,  -3),
                theta.t  = c(1, 1.1, 2.0,  -1.5),
                sigma    = sigma
                  )
kable(Power)

## ---- message=FALSE-----------------------------------------------------------
## Simulate data for treatment arm
data.trt <- rnorm(60, mean = 3, sd = sigma)

## first obtain posterior distributions...
post_SAM <- postmix(priormix = SAM.prior,   ## SAM Prior
                    data = data.crt)
post_trt <- postmix(priormix = unit_prior,  ## Non-informative prior
                    data = data.trt)

## Define the decision function
decision = decision2S(0.95, 0, lower.tail=FALSE)

## Decision-making
decision(post_trt, post_SAM)

## -----------------------------------------------------------------------------
sessionInfo()

## ----include=FALSE------------------------------------------------------------
options(.user_mc_options)

