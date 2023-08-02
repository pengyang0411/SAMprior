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

## ----results="asis",echo=FALSE------------------------------------------------
## Current trial data
ASAS20 <- data.frame(study = c('Baeten (2013)', 'Deodhar (2016)', 'Deodhar (2019)',
                               'Erdes (2019)', 'Huang (2019)', 'Kivitz (2018)',
                               'Pavelka (2017)', 'Sieper (2017)', 'Van der Heijde (2018)'),
                     n = c(6, 122, 104, 23, 153, 117, 76, 74, 87),
                     r = c(1, 35, 31, 10, 56, 55, 28, 21, 35))
kable(ASAS20)

## -----------------------------------------------------------------------------
# load R packages
library(ggplot2)
theme_set(theme_bw()) # sets up plotting theme
set.seed(22)
map_ASAS20 <- gMAP(cbind(r, n-r) ~ 1 | study,
                   family = binomial,
                   data = ASAS20, 
                   tau.dist = "HalfNormal", 
                   tau.prior = 1,
                   beta.prior = 2)

map_automix <- automixfit(map_ASAS20)
map_automix
plot(map_automix)$mix

## ---- message=FALSE-----------------------------------------------------------
n <- 35; r = 10 
wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.2,
                   n = n, r = r)
cat('SAM weight: ', wSAM)

## ---- message=FALSE-----------------------------------------------------------
wSAM <- SAM_weight(if.prior = map_automix, 
                   delta = 0.2,
                   method.w = 'PPR',
                   prior.odds = 3/7,
                   n = n, r = r)
cat('SAM weight: ', wSAM)

## ---- echo=FALSE, message=FALSE-----------------------------------------------
weight_grid <- seq(1, 34, by = 1)
weight_res  <- lapply(1:34, function(x) SAM_weight(if.prior = map_automix, 
                                                   delta = 0.2,
                                                   n = 35, r = x))
df_weight <- data.frame(grid   = weight_grid/35,
                        weight = unlist(weight_res))
qplot(grid, weight, data = df_weight, geom = "line", main= "SAM Weight") +
  xlab('Response Rate from Control trial')+ ylab('Weight') +
  geom_vline(xintercept = summary(map_automix)['mean'], linetype = 2, col = 'blue') 

## ---- message=FALSE-----------------------------------------------------------
SAM.prior <- SAM_prior(if.prior = map_automix, 
                       nf.prior = mixbeta(nf.prior = c(1,1,1)),
                       weight = wSAM)
SAM.prior

## ---- message=FALSE-----------------------------------------------------------
set.seed(123)
TypeI <- get_OC(if.prior = map_automix,       ## MAP prior from historical data
                nf.prior = mixbeta(c(1,1,1)), ## Non-informative prior for treatment arm
                delta    = 0.2,               ## CSD for SAM prior
                ## Method to determine the mixture weight for the SAM prior
                method.w = 'LRT',             
                n        = 35, n.t = 70,      ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,              ## Number of trials simulated
                if.MAP   = TRUE,              ## Output robust MAP prior for comparison
                weight   = 0.5,               ## Weight for robust MAP prior
                ## Response rates for control and treatment arms
                theta    = c(0.36, 0.36, 0.11, 0.55),
                theta.t  = c(0.34, 0.33, 0.11, 0.55)
                  )
kable(TypeI)

## ---- message=FALSE-----------------------------------------------------------
set.seed(123)
Power <- get_OC(if.prior = map_automix,       ## MAP prior from historical data
                nf.prior = mixbeta(c(1,1,1)), ## Non-informative prior for treatment arm
                delta    = 0.2,               ## CSD for SAM prior
                n        = 35, n.t = 70,      ## Sample size for control and treatment arms
                ## Decisions
                decision = decision2S(0.95, 0, lower.tail=FALSE), 
                ntrial   = 1000,              ## Number of trials simulated
                if.MAP   = TRUE,           ## Output robust MAP prior for comparison
                weight   = 0.5,               ## Weight for robust MAP prior
                ## Response rates for control and treatment arms
                theta    = c(0.37, 0.34, 0.16, 0.11),
                theta.t  = c(0.57, 0.54, 0.36, 0.31)
                  )
kable(Power)

## -----------------------------------------------------------------------------
## Sample size and number of responses for treatment arm
n_t <- 70; x_t <- 22 

## first obtain posterior distributions...
post_SAM <- postmix(priormix = SAM.prior,         ## SAM Prior
                    r = r,   n = n)
post_trt <- postmix(priormix = mixbeta(c(1,1,1)), ## Non-informative prior
                    r = x_t, n = n_t)

## Define the decision function
decision = decision2S(0.95, 0, lower.tail=FALSE)

## Decision-making
decision(post_trt, post_SAM)

## -----------------------------------------------------------------------------
sessionInfo()

## ----include=FALSE------------------------------------------------------------
options(.user_mc_options)

