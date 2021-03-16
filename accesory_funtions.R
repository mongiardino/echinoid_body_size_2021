`%not in%` = function(x, table) is.na(match(x, table, nomatch=NA_integer_))

#equation to calculate AICc from AIC
aic_to_aicc = function(aic = NULL, param = NULL, tree = NULL) {
  return(aic + (2*param*(param+1)/(length(tree$tip.label)-param-1)))
}

##Pulsed evolution from Schraiber & Landis (2018)-------------------------------------------------------
#Everything from here until the end of the file 

#some packages that this needs

packages <- c('fBasics', 'neldermead', 'gsl', 'stabledist', 'statmod')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

require(fBasics) #to simulate alpha stable guys, ?dstable
require(neldermead)
require(gsl) #for better bessel functions
require(stabledist)
require(statmod)

##########some characteristic exponents#########

get_cf = function(m) {
  if (m=="BM")    return(ebm)
  if (m=="OU")    return(ebm)
  if (m=="EB")    return(ebm)
  if (m=="JN")    return(ejn)
  if (m=="VG")    return(evg)
  if (m=="NIG")   return(enig)
  if (m=="BMJN")  return(ebmjn)
  if (m=="BMVG")  return(ebmvg)
  if (m=="BMNIG") return(ebmnig)
  if (m=="EBJN")  return(ejn)
  if (m=="EBVG")  return(evg)
  if (m=="EBNIG") return(enig)
}


####################
# GAUSSIAN PROCESS #
####################

# BM: Brownian motion
ebm = function(k,t,sigma_bm,singularity=FALSE,variance=FALSE) {
  if (singularity) {
    return(c(-Inf,Inf))
  }
  if (variance) {
    return(sigma_bm^2)
  }
  if (sigma_bm < 0) {
    -Inf	
  } else {
    -sigma_bm^2*k^2/2*t
  }
}


####################
# COMPOUND POISSON #
####################

# JN: compound Poisson process w/ normal jumps
ejn = function(k,t,lambda_jn,sigma_jn,singularity=FALSE,variance=FALSE) {
  if (singularity) {
    return(c(-Inf,Inf))
  }
  if (variance) {
    return(lambda_jn*sigma_jn^2)
  }
  if (lambda_jn < 0 || sigma_jn < 0){
    -Inf	
  } else {
    lambda_jn*t*(exp(-sigma_jn^2*k^2/2)-1)
  }
}

# BM+JN: Brownian motion + compound Poisson process w/ normal jumps
ebmjn = function(k,t,sigma_bm,lambda_jn,sigma_jn,singularity=FALSE,variance=FALSE) {
  if (singularity) {
    return(c(-Inf,Inf))
  }
  if (variance) {
    return(sigma_bm^2+lambda_jn*sigma_jn^2)
  }
  if (sigma_bm < 0 || lambda_jn < 0 || sigma_jn < 0){
    -Inf	
  } else {
    -sigma_bm^2*k^2/2*t + lambda_jn*t*(exp(-sigma_jn^2*k^2/2)-1)
  }
}


##################
# VARIANCE GAMMA #
##################

# VG: variance Gamma process
evg = function(k,t,sigma_vg,nu_vg,mu_vg=0,singularity=FALSE,variance=FALSE) {
  #parameterized as in Madan, Carr and Chang (1998)
  # sigma_vg : \sigma, volatility of BM
  # nu_vg    : \nu, variance rate of gamma time change
  # mu_vg    : \theta, directional trend
  # k        : u, the cf variable
  if (variance) {
    return(sigma_vg^2)
  }
  if (singularity) {
    v = sqrt(as.complex(-mu_vg^2-2*sigma_vg^2/nu_vg))
    return(c(1i*mu_vg+v,1i*mu_vg-v)/sigma_vg^2)
  }
  if (nu_vg < 0 || sigma_vg < 0){
    -Inf
  } else {
    -t/nu_vg*log(1-(1i)*mu_vg*nu_vg*k + k^2*sigma_vg^2*nu_vg/2 )
  }
}

# BM+VG: Brownian motion + variance Gamma process
ebmvg = function(k,t,sigma_bm,sigma_vg,nu_vg,mu_vg=0,singularity=FALSE,variance=FALSE) {
  #parameterized as in Madan, Carr and Chang (1998)
  # sigma_vg : \sigma, volatility of BM
  # nu_vg    : \nu, variance rate of gamma time change
  # mu_vg    : \theta, directional trend
  # k        : u, the cf variable
  if (variance) {
    return(sigma_bm^2+sigma_vg^2)
  }
  if (singularity) {
    v = sqrt(as.complex(-mu_vg^2-2*sigma_vg^2/nu_vg))
    return(c(1i*mu_vg+v,1i*mu_vg-v)/sigma_vg^2)
  }
  if (sigma_bm < 0 || nu_vg < 0 || sigma_vg < 0){
    -Inf
  } else {
    -sigma_bm^2*k^2/2*t - t/nu_vg*log(1-(1i)*mu_vg*nu_vg*k + k^2*sigma_vg^2*nu_vg/2 )
  }
}


###########################
# NORMAL INVERSE GAUSSIAN #
###########################

# NIG: normal inverse Gaussian process
enig = function(k,t,alpha_nig,delta_nig,beta_nig=0,mu_nig=0,singularity=FALSE,variance=FALSE) {
  if (variance) {
    return(t*alpha_nig^2*delta_nig/(alpha_nig^2-beta_nig^2)^(3/2))
  }
  if (singularity) {
    return(c(-Inf,Inf))
  }
  
  if (alpha_nig < abs(beta_nig)) {
    -Inf
  } else if (delta_nig < 0) {
    -Inf
  } else { 
    gamma_nig = sqrt(alpha_nig^2-beta_nig^2)
    t*((1i)*mu_nig*k+delta_nig*(gamma_nig-sqrt(alpha_nig^2-(beta_nig+(1i)*k)^2)))
  }
}

# BM+NIG: Brownian motion + normal inverse Gaussian process
ebmnig = function(k,t,sigma_bm,alpha_nig,delta_nig,beta_nig=0,mu_nig=0,singularity=FALSE,variance=FALSE) {
  if (variance) {
    return(sigma_bm^2 + t*alpha_nig^2*delta_nig/(alpha_nig^2-beta_nig^2)^(3/2))
  }
  if (singularity) {
    return(c(-Inf,Inf))
  }
  if (sigma_bm < 0) {
    -Inf
  } else if (alpha_nig < abs(beta_nig)) {
    -Inf
  } else if (delta_nig < 0) {
    -Inf
  } else { 
    gamma_nig = sqrt(alpha_nig^2-beta_nig^2)
    -sigma_bm^2*k^2/2*t + t*((1i)*mu_nig*k+delta_nig*(gamma_nig-sqrt(alpha_nig^2-(beta_nig+(1i)*k)^2)))
  }
}

#source("levy_pruning_prob.r")
#source("levy_pruning_cf.r")

################
# OPTIMIZATION #
################

# ensures initial values are within tolerance limits
tol_safe = function(x, tol) {
  # some redundancy, but modify if needed
  if (tol < 0 && x > tol) {
    x = tol - rexp(1,100)
  } else if (tol > 0 && x < tol) {
    x = tol + rexp(1,100)
  }
  return(x)
}

# Arguments
# phy               : an ape tree object
# dat               : a vector of univariate trait data with names
#                     matching phy$tip.label
# model             : a string indicating the model (BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG)
# par               : initial parameters (must match model dimensions/bounds)
# sigma_tip         : estimate sigma_tip value? (default==TRUE)
# tol               : tolerance to terminate optimization (when delta-lnL or delta-param < tol)
# maxfun            : max number of function calls during optimization
# maxiter           : max number of iterations during optimization
# weight_scaling    : weight scaling to improve integration
# silent            : do not print text during optimization
#
# Returns list object with the following elements
# params            : MLE parameters
# n_params          : number of parameters
# lnL               : MLE log likelihood
# AIC               : Akiake Information Criterion score
# optim             : neldermead optimization object

fit_reml_levy = function(phy,dat,model,par=NA,sigma_tip=T,tol=1e-4,maxfun=2000,maxiter=1000,weight_scaling=1e5,silent=T) {
  
  # validate model choice
  model = toupper(model)
  if ( !(model %in% c("BM","OU","EB","JN","VG","NIG","BMJN","BMVG","BMNIG","EBJN","EBVG","EBNIG")) ) {
    stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG, EBJN, EBVG, EBNIG")
  }
  
  # model tip noise?
  tip_noise = NA # estimate tipnoise
  tip_noise_lower = 1e-9
  if (!sigma_tip) {
    tip_noise = 0
    tip_noise_lower = 0
  }
  
  # initial parameter guess
  s = exp( runif(3, -0.5, 0.5) )
  alpha_ou = 1e-4
  decay_eb = -1e-4
  if (model=="OU" && runif(1) < 0.5) {
    # OU
    phy2 = chronos(phy, quiet=T)
    par_init = fitContinuous(phy2,dat,model="OU",SE=tip_noise)$opt
    alpha_ou = tol_safe( par_init$alpha * rbeta(1,1,5), tol )
  } else if (model%in%c("EB","EBJN","EBVG","EBNIG") && runif(1) < 0.5) {
    # EB
    phy2 = chronos(phy, quiet=T)
    par_init = fitContinuous(phy2,dat,model="EB",SE=tip_noise)$opt
    decay_eb = tol_safe( par_init$a * rbeta(1,1,5), -tol )
  } else {
    # all other models
    par_init = fitContinuous(phy,dat,model="BM",SE=tip_noise)$opt
  }
  sig_bm = tol_safe( sqrt(par_init$sigsq) * s[1], tol )
  if (sigma_tip) {
    sig_tip = tol_safe( par_init$SE * s[2], tol )
  } else {
    sig_tip = 0
  }
  proc_var = sig_bm^2
  if (runif(1) < 0.5) {
    proc_kurt = rgamma(1, shape=2, rate=0.2)
  } else {
    proc_kurt = runif(1,0.5,5)
  }
  frac_of_var = rbeta(1,2,2)
  
  
  # BROWNIAN MOTION
  if (model=="BM") {
    # set initial params
    if (any(is.na(par))) {
      par = c(sig_bm, sig_tip)
    }
    lower = rep(1e-9,2)
    upper = rep(1e+3,2)
    lower[2] = tip_noise_lower
    
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[2] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"BM",sigma_bm=x[1],theta_ou=0,decay_eb=0,sigma_tip=x[2],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # ORNSTEIN-UHLENBECK
  #else if (all.equal(phi,ebm)==T && is.ou && !is.eb) {
  else if (model=="OU") {
    # set initial params
    if (any(is.na(par))) {
      par = c(sig_bm, alpha_ou, sig_tip)
    }
    lower = rep(1e-9,3)
    upper = rep(1e+3,3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[3] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"OU",sigma_bm=x[1],theta_ou=x[2],decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # EARLY BURST
  else if (model=="EB") {
    # set initial params
    if (any(is.na(par))) {
      par = c(sig_bm, decay_eb, sig_tip)
    }
    lower = c(1e-9,-1e+3,1e-9)
    upper = c(1e+3,-1e-9,1e+3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[3] = 0
      if (x[1] < 0 || x[2] > 0 || x[3] < 0) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"EB",sigma_bm=x[1],theta_ou=0,decay_eb=x[2],sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # BROWNIAN MOTION + COMPOUND POISSON WITH NORMAL JUMPS
  else if (model=="BMJN") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_JN(proc_var, proc_kurt, frac_of_var)
      sigma_bm  = par_tmp$sigma.bm
      lambda_jn = par_tmp$lambda.jn
      delta_jn  = par_tmp$delta.jn
      par = c(sigma_bm, lambda_jn, delta_jn, sig_tip)
    }
    lower = rep(1e-9,4)
    upper = rep(1e+6,4)
    lower[4] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"BMJN",sigma_bm=x[1],lambda_jn=x[2],sigma_jn=x[3],
                             theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # PURE COMPOUND POISSON WITH NORMAL JUMPS
  else if (model=="JN") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_JN(proc_var, proc_kurt, 1)
      lambda_jn = par_tmp$lambda.jn
      delta_jn  = par_tmp$delta.jn
      par = c(lambda_jn, delta_jn, sig_tip)
    }
    lower = rep(1e-9,3)
    upper = rep(1e+3,3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[3] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"JN",lambda_jn=x[1],sigma_jn=x[2],
                             theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # BROWNIAN MOTION + VARIANCE GAMMA
  else if (model=="BMVG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_VG(proc_var, proc_kurt, frac_of_var)
      par = c(par_tmp$sigma.bm, par_tmp$nu.vg, par_tmp$sigma.vg, sig_tip)
    }
    lower = rep(1e-9,4)
    upper = rep(1e+3,4)
    lower[4] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"BMVG",sigma_bm=x[1],sigma_vg=x[2],nu_vg=x[3],
                             theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # PURE VARIANCE GAMMA
  else if (model=="VG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_VG(proc_var, proc_kurt, 1)
      par = c(par_tmp$nu.vg, par_tmp$sigma.vg, sig_tip)
    }
    lower = rep(1e-9,3)
    upper = rep(1e+3,3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[3] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"VG",sigma_vg=x[1],nu_vg=x[2],
                             theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # BROWNIAN MOTION + NORMAL INVERSE GAUSSIAN
  else if (model=="BMNIG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_NIG(proc_var, proc_kurt, frac_of_var)
      par = c(par_tmp$sigma.bm, par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
    }
    lower = rep(1e-9,4)
    upper = rep(1e+3,4)
    lower[4] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"BMNIG",sigma_bm=x[1],alpha_nig=x[2],delta_nig=x[3],
                             theta_ou=0,decay_eb=0,sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # PURE NORMAL INVERSE GAUSSIAN
  else if (model=="NIG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_NIG(proc_var, proc_kurt, frac_of_var)
      par = c(par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
    }
    lower = rep(1e-9,3)
    upper = rep(1e+3,3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[3] = 0
      if (any(x<0)) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"NIG",alpha_nig=x[1],delta_nig=x[2],
                             theta_ou=0,decay_eb=0,sigma_tip=x[3],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  # EARLY BURST w/ RATE DECAYING JN
  else if (model=="EBJN") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_JN(proc_var, proc_kurt, 1)
      lambda_jn = par_tmp$lambda.jn
      delta_jn  = par_tmp$delta.jn
      par = c(decay_eb, lambda_jn, delta_jn, sig_tip)
    }
    lower = c(-1e+3,1e-9,1e-9,1e-9)
    upper = c(-1e-9,1e+3,1e+3,1e+3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"EBJN",theta_ou=0,decay_eb=x[1],lambda_jn=x[2],sigma_jn=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  
  # EARLY BURST w/ TIME DECAYING VG
  else if (model=="EBVG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_VG(proc_var, proc_kurt, 1)
      par = c(decay_eb, par_tmp$nu.vg, par_tmp$sigma.vg, 1)
    }
    lower = c(-1e+3,1e-9,1e-9,1e-9)
    upper = c(-1e-9,1e+3,1e+3,1e+3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"EBVG",theta_ou=0,decay_eb=x[1],nu_vg=x[2],sigma_vg=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  
  # EARLY BURST w/ TIME DECAYING NIG
  else if (model=="EBNIG") {
    # set initial params
    if (any(is.na(par))) {
      par_tmp = get_params_NIG(proc_var, proc_kurt, 1)
      par = c(decay_eb, par_tmp$alpha.nig, par_tmp$delta.nig, sig_tip)
    }
    lower = c(-1e+3,1e-9,1e-9,1e-9)
    upper = c(-1e-9,1e+3,1e+3,1e+3)
    lower[3] = tip_noise_lower
    fn = function(x=NULL, index=NULL, fmsfundata=NULL) {
      if (!sigma_tip) x[4] = 0
      if (x[1] > 0 || x[2] <  0 || x[3] < 0 || x[4] < 0) {
        ret = -Inf
      } else {
        ret = llik_reml_levy(phy,dat,"EBNIG",theta_ou=0,decay_eb=x[1],alpha_nig=x[2],delta_nig=x[3],sigma_tip=x[4],weight_scaling=weight_scaling,silent=silent)
      }
      if (any(ret == -Inf)) ret = -(0.1/tol)
      return(list(f=-sum(ret),
                  g=c(),
                  c=c(),
                  gc=c(),
                  index=index,
                  this=list(costfargument=fmsfundata)))
    }
  }
  else {
    return(NA)
  }
  
  x0 = transpose(par)
  nm = neldermead()
  nm = neldermead.set(nm,'numberofvariables',length(par))
  nm = neldermead.set(nm,'function',fn)
  nm = neldermead.set(nm,'x0',x0)
  nm = neldermead.set(nm,'verbose',FALSE)
  nm = neldermead.set(nm,'storehistory',TRUE)
  nm = neldermead.set(nm,'verbosetermination',FALSE)
  nm = neldermead.set(nm,'method','box')
  nm = neldermead.set(nm,'boundsmin',lower)
  nm = neldermead.set(nm,'boundsmax',upper)
  nm = neldermead.set(nm,'tolfunmethod', TRUE)
  nm = neldermead.set(nm,'tolfunrelative',tol)
  nm = neldermead.set(nm,'tolxmethod',TRUE)
  nm = neldermead.set(nm,'tolxrelative',tol)
  nm = neldermead.set(nm,'maxfunevals',maxfun)
  nm = neldermead.set(nm,'maxiter',maxiter)
  nm = neldermead.set(nm,'boxtermination',TRUE)
  nm = neldermead.set(nm,'boxtolf',tol*1e-2)
  nm = neldermead.search(nm)
  nm$gc = gc()
  nm$phy = phy
  
  if (!sigma_tip) {
    tip_noise_idx = length( nm$optbase$xopt[,1] )
    p = nm$optbase$xopt[,1]
    p[length(p)]=0
    nm$optbase$xopt[,1] = p
  }
  
  if (!silent) cat ("...done!\n")
  
  # return
  results = list()
  results$model = model
  results$dat = dat
  results$phy = phy
  results$params = format_params(nm$optbase$xopt[,1], model)
  results$n_params = sum( results$params != 0 )
  results$lnL = -nm$optbase$fopt
  results$AIC = 2*(results$n_params-results$lnL)
  results$optim = nm
  return(results)
}


get_best_idx = function(x) {
  n = length(x)
  best = which.min(unlist(sapply(1:n,
                                 function(i)
                                 {
                                   z=x[[i]]$optbase$fopt
                                   if (!is.numeric(z))
                                     z = NA
                                   return(z)
                                 }
  )))
  return(best)
}

#source("levy_pruning_tools.r")

rescale_error_cf = function(left_scale,right_scale,left_times,right_times,sigma_left,sigma_right,t_l,t_r,v_l,v_r,sig_p2) {
  #rescale left by multiplying by RIGHT branch rescaling factor
  new_left_scale = (sig_p2*t_r+v_r)*c(1,left_scale)/
    (sig_p2*(t_l+t_r)+v_l+v_r)
  #rescale right by multiplying by LEFT branch rescaling factor
  new_right_scale = (sig_p2*t_l+v_l)*c(1,right_scale)/
    (sig_p2*(t_l+t_r)+v_l+v_r)
  #add the new branch lengths
  new_left_times = c(t_l,left_times)
  new_right_times = c(t_r,right_times)
  #add the new sigmas
  new_sigma_left = c(0,sigma_left)
  new_sigma_right = c(0,sigma_right)
  return(list(scales=c(new_left_scale,new_right_scale),times=c(new_left_times,new_right_times),sigma=c(new_sigma_left,new_sigma_right)))
}

#NB: this function expects phi to be the characteristic exponent!
error_cf.vectorized = function(k,scale,times,phi,...) {
  if (length(times)!=length(scale)) {
    stop("length(times) = ",length(times)," != length(scale) = ",length(scale))
  } 
  if (length(times) == 0) {
    return(rep(1,length(k)))
  }
  #vectorize the computation
  scaleKmat = scale%*%t(k)
  timeMat = matrix(nrow=length(times),ncol=length(k),times)
  errMat = phi(scaleKmat,timeMat,...)
  cur_err = colSums(errMat)
  return(cur_err)
}

#NB: This function expects phi to be the characteristic exponent!
error_cf = function(k,sigma_tip,scale,times,phi,...) {
  
  # here be the tipward pass to compute the cf
  cur_err = 0
  #TODO: THIS FLAGGING ISN'T WORKING
  for (i in 1:length(scale)) {
    #TODO: for the tips, need to add a *different* phi
    if (sigma_tip[i] == 0) {
      cur_err = cur_err + phi(scale[i]*k,times[i],...)
    } else {
      #cur_err = cur_err - sigma_tip[i]^2*scale[i]^2*k^2/2
      cur_err = cur_err - sigma_tip[i]*sigma_tip[i]*scale[i]*scale[i]*k*k/2
    }
  }
  return(cur_err)
}

#NB: This function expects phi to be the characteristic exponent!
contrast_cf = function(k,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...) {
  #compute the cf of the change to the left
  cur_cf = phi(k,t_l,...)
  #compute the cf of the change to the right (evaluated at negative)
  cur_cf = cur_cf+phi(-k,t_r,...)
  #compute the cf of the error to the left
  cur_cf = cur_cf+error_cf(k,sigma_left,left_scale,left_times,phi,...)
  #compute the cf of the error to the right (evaluated at negative)
  cur_cf = cur_cf+error_cf(-k,sigma_right,right_scale,right_times,phi,...)
  return(exp(cur_cf))
}

weight = function(k,p=8,q=4) {
  # returns 0.5*erfc(k/p-q)
  return(pnorm(-sqrt(2)*(k/p-q)))
}

fourier_real = function(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,..., weight_scaling=100000) {
  #gets the real part of the fourier integral
  return(function(k) {
    cur_cf = contrast_cf(k,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...)
    return(weight(k,sqrt(weight_scaling/abs(x)),sqrt(weight_scaling*abs(x)/4))*(Re(cur_cf)*cos(-k*x)-Im(cur_cf)*sin(-k*x)))
  })
}

contrast_likelihood = function(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...,weight_scaling=100000) {
  contrast_re = fourier_real(x,sigma_left,sigma_right,left_scale,right_scale,left_times,right_times,t_l,t_r,phi,...,weight_scaling=weight_scaling)
  #compute integral
  real_part = 1/(pi)*integrate(contrast_re,0,Inf,subdivisions=1000000L,rel.tol=1e-50,abs.tol=1e-50,stop.on.error=F)$value
  return(real_part)
}

llik_reml_levy = function(phy,dat,model,...,return_inf=TRUE,theta_ou=0,decay_eb=0,sigma_tip=0,weight_scaling=100000,silent=T) {
  
  ### INITIALIZE DATA ###
  #sort the tree so that it can be pruned
  phy = reorder(phy,"postorder")
  
  #ensure that data is a named vector
  if (!is.vector(dat)) {
    stop("ERROR: dat is not a vector")
  }
  if (is.null(names(dat))) {
    stop("ERROR: dat does not have names")
  } 
  # label and order data
  old_length = length(dat)
  dat = dat[phy$tip.label]
  if (length(dat)!=old_length) {
    warning("There are discrepancies between dat and phy")
  }
  
  ### INITIALIZE MODEL ###
  model = toupper(model)
  if ( !(model %in% c("BM","OU","EB","JN","VG","NIG","BMJN","BMVG","BMNIG","EBJN","EBVG","EBNIG")) ) {
    stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, BMJN, BMVG, BMNIG, EBJN, EBVG, EBNIG")
  }
  phi = get_cf(model)
  
  # rescale by EB/OU parameters for use with ebm
  if (theta_ou > 0 && model=="OU") {
    phy = OU.brlen(phy=phy, theta=theta_ou)
  }
  else if (decay_eb != 0 && model%in%c("EB","EBJN","EBVG","EBNIG")) {
    phy = EB.brlen(phy=phy, r=decay_eb)
  }
  
  ### INITIALIZE WORKSPACE ####
  scales = list()
  times = list()
  sigma_tip_list = list()
  pseudo_obs = rep(0,Nnode(phy,internal.only=FALSE))
  contrast_lnL = rep(0,Nnode(phy,internal.only=FALSE))
  contrast_val = rep(0,Nnode(phy,internal.only=FALSE))
  error_var = rep(0,Nnode(phy,internal.only=FALSE))
  #get the process variance
  sig_p2 = phi(1,1,...,variance=TRUE)
  
  
  ### COMPUTE LIKELIHOOD ###
  #first initialize the tips
  for (i in 1:length(phy$tip.label)) {
    scales[[i]] = 1 #initalize this at 1 so that the tip errors are included
    times[[i]] = 0
    sigma_tip_list[[i]] = sigma_tip
    pseudo_obs[i] = dat[i]
    contrast_lnL[i] = 0
    error_var[i] = sigma_tip^2
  }
  #then loop over the internal nodes
  for (i in seq(1,nrow(phy$edge),2)) {
    #get the nodes
    cur_node = phy$edge[i,1]
    cur_left = phy$edge[i,2]
    cur_right = phy$edge[i+1,2]
    #get the relevant scalings and times
    left_scale = scales[[cur_left]]
    right_scale = scales[[cur_right]]
    sigma_left = sigma_tip_list[[cur_left]]
    sigma_right = sigma_tip_list[[cur_right]]
    left_time = times[[cur_left]]
    right_time = times[[cur_right]]
    t_l = phy$edge.length[i]
    t_r = phy$edge.length[i+1]
    v_l = error_var[cur_left]
    v_r = error_var[cur_right]
    #compute the contrast likelihood
    contrast = pseudo_obs[cur_left]-pseudo_obs[cur_right]
    cur_like = contrast_likelihood(contrast,sigma_left,sigma_right,left_scale,right_scale,left_time,right_time,t_l,t_r,phi,...,weight_scaling=weight_scaling)
    contrast_val[cur_node] = contrast
    if (cur_like > 0) {
      contrast_lnL[cur_node] = log(cur_like)			
    } else if (!return_inf) {
      contrast_lnL[cur_node] = -300
    } else {
      return(-Inf)
    }
    #compute the new error var
    new_error_var = (sig_p2*t_l+v_l)*(sig_p2*t_r+v_r)/(sig_p2*(t_l+t_r)+v_l+v_r)
    error_var[cur_node] = new_error_var
    #compute the new error scaling
    new_error = rescale_error_cf(left_scale,right_scale,left_time,right_time,sigma_left,sigma_right,t_l,t_r,v_l,v_r,sig_p2)
    scales[[cur_node]] = new_error$scales
    times[[cur_node]] = new_error$times
    sigma_tip_list[[cur_node]] = new_error$sigma
    #compute the new pseudo-observation
    pseudo_obs[cur_node] = ((sig_p2*t_r+v_r)*pseudo_obs[cur_left]+
                              (sig_p2*t_l+v_l)*pseudo_obs[cur_right])/
      (sig_p2*(t_l+t_r)+v_l+v_r)
    #pseudo_obs[cur_node] = pseudo_obs[cur_node]/(t_l+t_r)
  }
  
  results = c(sum(contrast_lnL),...)
  names(results)[1] = "lnL"
  if (decay_eb!=0) {
    names(decay_eb)="decay_eb"
    results = c(results, decay_eb)
  }
  if (theta_ou>0) {
    names(theta_ou)="theta_ou"
    results = c(results, theta_ou)
  }
  names(sigma_tip)="sigma_tip"
  results = c(results,sigma_tip)
  
  if (!silent) print(results)
  
  return(contrast_lnL)
}

#source("levy_pruning_tools.r")

rlevy = function(phy, model, par, n=1) {
  
  # format data
  par = check_args(model, par)
  print(par)
  phy = reorder(phy, "cladewise")
  
  # simulate
  if (model=="BM") {
    x = sim_bm(phy=phy,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
  } else if (model=="OU") {
    phy2 = OU.brlen(phy=phy,theta=par$theta_ou)
    phy2 = reorder(phy2, "cladewise")
    x = sim_bm(phy=phy2,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
  } else if (model=="EB") {
    phy2 = EB.brlen(phy=phy, r=par$decay_eb)
    phy2 = reorder(phy2, "cladewise")
    x = sim_bm(phy=phy2,sigma_bm=par$sigma_bm,sigma_tip=par$sigma_tip)
  } else if (model=="JN") {
    x = sim_cpp(phy=phy,
                sigma=0,
                lambda=par$lambda_jn,
                kernel=rnorm,
                mean=0,
                sd=par$delta_jn,
                sigma_tip=par$sigma_tip)
  } else if (model=="VG") {
    x = sim_vg(phy=phy,
               sigma_bm=0,
               sigma_vg=par$sigma_vg,
               nu_vg=par$nu_vg,
               mu_vg=0,
               sigma_tip=par$sigma_tip)
  } else if (model=="NIG") {
    x = sim_nig(phy=phy,
                sigma_bm=0,
                alpha_nig=par$alpha_nig,
                beta_nig=0,
                delta_nig=par$delta_nig,
                mu_nig=0,
                sigma_tip=par$sigma_tip) 
  } else if (model=="BMJN") {
    x = sim_cpp(phy=phy,
                sigma=par$sigma_bm,
                lambda=par$lambda_jn,
                kernel=rnorm,
                mean=0,
                sd=par$delta_jn,
                sigma_tip=par$sigma_tip)
  } else if (model=="BMVG") {
    x = sim_vg(phy=phy,
               sigma_bm=par$sigma_bm,
               sigma_vg=par$sigma_vg,
               nu_vg=par$nu_vg,
               mu_vg=0,
               sigma_tip=par$sigma_tip)
  } else if (model=="BMNIG") {
    x = sim_nig(phy=phy,
                sigma_bm=par$sigma_bm,
                alpha_nig=par$alpha_nig,
                beta_nig=0,
                delta_nig=par$delta_nig,
                mu_nig=0,
                sigma_tip=par$sigma_tip)
  } else if (model=="AS") {
    x = sim_stable(phy=phy,
                   sigma=0,
                   alpha=par$alpha_as,
                   cee=par$c_as,
                   sigma_tip=par$sigma_tip)
  } else if (model=="BMAS") {
    x = sim_stable(phy=phy,
                   sigma=par$sigma_bm,
                   alpha=par$alpha_as,
                   cee=par$c_as,
                   sigma_tip=par$sigma_tip)
  }
  
  return(x)
}


jump_dif = function(sigma,lambda,kernel,...) {
  #sigma = rate*t of brownian motion
  #lambda = rate*t of poisson process
  #kernel = function to draw from jump kernel
  #... = parameters for kernel (MUST BE IN CORRECT ORDER!)
  numJumps = rpois(1,lambda)
  #Jbranch <<- append(Jbranch,numJumps)
  curState = 0
  curState = rnorm(1,0,sd=sigma)
  if (numJumps > 0) {
    curState = curState + sum(kernel(numJumps,...))
  }
  return(curState)	
}

sim_bm = function(phy, sigma_bm, sigma_tip = 0) {
  #phy is an ape-format tree
  #sigma is rate of brownian motion
  #nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
  nodes = rep(NA, nrow(phy$edge)+1)
  nodes[phy$edge[1,1]] = 0 # root value
  for (i in 1:length(phy$edge[,1])) {
    curLen = phy$edge.length[i]
    dx = rnorm(n=1, mean=0, sd=sigma_bm*sqrt(curLen)) #,lambda*curLen,kernel,...)
    nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + dx
    names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
  }
  nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
  return( nodes + rnorm(length(phy$tip.label),mean=0,sd=sigma_tip) )
}

sim_stable = function(phy, sigma, alpha, cee, sigma_tip=0) {
  #phy is ape-format tree
  #sigma is rate of brownian motion
  #cee is the scale parameter
  #alpha is the stability parameter
  nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
  nodes[phy$edge[1,1]] = 0
  for (i in 1:length(phy$edge[,1])) {
    curLen = phy$edge.length[i]
    nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, mean=0, sd = sigma*sqrt(curLen)) + rstable(1,alpha=alpha,beta=0,gamma=(curLen)^(1/alpha)*cee)
    names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
  }
  #print(nodes)
  nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
  return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
  
}

sim_vg = function(phy,sigma_bm,sigma_vg,nu_vg,mu_vg,sigma_tip=0) {
  #parameterized as in Madan, Carr, Chang (1998)
  #phy is ape-format tree
  #sigma_bm is brownian motion rate
  #nu_vg controls kurtosis of VG
  #sigma_vg controls rate of VG
  #mu_vg controls skewness of VG
  
  #implements as a time-changed brownian motion.
  nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
  nodes[phy$edge[1,1]] = 0
  for (i in 1:length(phy$edge[,1])) {
    curLen = phy$edge.length[i]
    curTimeChange = rgamma(1,curLen/nu_vg,scale=nu_vg)
    nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma_bm*sqrt(curLen)) + rnorm(1, sd = sigma_vg*sqrt(curTimeChange), mean = mu_vg*curTimeChange)
    names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
  }
  #return(rnorm(length(phy$tip.label),nodes[1:length(phy$tip.label)],sigma_tip))
  nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
  return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

sim_nig = function(phy,sigma_bm,alpha_nig,beta_nig,delta_nig,mu_nig,sigma_tip=0) {
  #parametrized as in PhD seminar for modeling normal inverse gaussian processses
  #www2.math.uni-wuppertal.de/~ruediger/pages/vortraege/ws1112/nig_5.pdf
  #phy is ape-format tree
  #sigma_bm is brownian motion rate
  #alpha_nig controls kurtosis of NIG
  #beta_nig controls skewness of NIG
  #delta_nig controls the scale of NIG
  #mu_nig is location parameter of NIG
  
  #implements as time-changed brownian motion
  gamma_nig = sqrt(alpha_nig^2 - beta_nig^2)
  nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
  nodes[phy$edge[1,1]] = 0
  for (i in 1:length(phy$edge[,1])) {
    curLen = phy$edge.length[i]
    curTimeChange = rinvgauss(1,delta_nig/gamma_nig*curLen,delta_nig^2*curLen^2)
    nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + rnorm(1, sd = sigma_bm*sqrt(curLen)) + rnorm(1, sd = sqrt(curTimeChange), mean = mu_nig*curLen + beta_nig*curTimeChange)
    names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
  }
  nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
  return(nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

rNIG = function(n=1,sigma_bm,alpha_nig,beta_nig,delta_nig,mu_nig,t=1) {
  #parametrized as in PhD seminar for modeling normal inverse gaussian processses
  #www2.math.uni-wuppertal.de/~ruediger/pages/vortraege/ws1112/nig_5.pdf
  #phy is ape-format tree
  #sigma_bm is brownian motion rate
  #alpha_nig controls kurtosis of NIG
  #beta_nig controls skewness of NIG
  #delta_nig controls the scale of NIG
  #mu_nig is location parameter of NIG
  
  #implements as time-changed brownian motion
  gamma_nig = sqrt(alpha_nig^2 - beta_nig^2)
  curLen = t
  curTimeChange = rinvgauss(n,delta_nig/gamma_nig*curLen,delta_nig^2*curLen^2)
  x = rnorm(n, sd=sigma_bm*sqrt(curLen)) + rnorm(n, sd=sqrt(curTimeChange), mean=mu_nig*curLen + beta_nig*curTimeChange)
  return(x)
}

#this simulates a compound poisson model with jump kernel given by kernel and parameters given by ...
sim_cpp = function(phy, sigma, lambda, kernel,...,sigma_tip = 0) {
  #phy is an ape-format tree
  #sigma is rate of brownian motion
  #lambda is rate of poisson process
  #... are kernel parameters
  #nodes = vector(length=(phy$Nnode+length(phy$tip.label)))
  nodes = rep(NA, nrow(phy$edge)+1)
  nodes[phy$edge[1,1]] = 0 #runif(1,-5,5) #sets the root
  for (i in 1:length(phy$edge[,1])) {
    curLen = phy$edge.length[i]
    dx = jump_dif(sigma*sqrt(curLen),lambda*curLen,kernel,...)
    nodes[phy$edge[i,2]] = nodes[phy$edge[i,1]] + dx
    names(nodes)[phy$edge[i,2]] = phy$tip.label[phy$edge[i,2]]
    
  }
  nodes = (nodes[1:length(phy$tip.label)])[phy$tip.label]
  return( nodes + rnorm(length(phy$tip.label),0,sigma_tip) )
}

check_args = function(model, par) {
  model = toupper(model)
  if ( !(model %in% c("BM","OU","EB","JN","VG","AS","NIG","BMJN","BMVG","BMNIG","BMAS")) ) {
    stop("Provided model type invalid. Valid models: BM, OU, EB, JN, VG, NIG, AS, BMJN, BMVG, BMNIG, BMAS")
  }
  
  if (!("sigma_tip" %in% names(par)))
    par$sigma_tip = 0
  
  # things with a Brownian component
  if (model=="BM" || model=="OU" || model=="EB" || 
      model=="BMJN" || model=="BMVG" || model=="BMNIG" || model=="BMAS") {
    
    if (!("sigma_bm" %in% names(par)))
      stop("Value for par$sigma_bm missing!")
  }
  
  # Ornstein-Uhlenbeck model
  if (model=="OU") {
    if (!("theta_ou" %in% names(par)))
      stop("Value for par$theta_ou missing!")
  }
  
  # Early Burst model
  if (model=="EB") {
    if (!("decay_eb" %in% names(par)))
      stop("Value for par$decay_eb missing!")
  }
  
  # jump normal models
  if (model=="JN" || model=="BMJN") {
    if (!("lambda_jn" %in% names(par)))
      stop("Value for par$lambda_jn missing!")
    if (!("delta_jn" %in% names(par)))
      stop("Value for par$delta_jn missing!")
  }
  
  # variance gamma models
  if (model=="VG" || model=="BMVG") {
    if (!("nu_vg" %in% names(par)))
      stop("Value for par$nu_vg missing!")
    if (!("sigma_vg" %in% names(par)))
      stop("Value for par$sigma_vg missing!")
  }
  
  # normal inverse gaussian models
  if (model=="NIG" || model=="BMNIG") {
    if (!("alpha_nig" %in% names(par)))
      stop("Value for par$alpha_nig missing!")
    if (!("delta_nig" %in% names(par)))
      stop("Value for par$delta_nig missing!")
  }
  
  # alpha stable models
  if (model=="AS" || model=="BMAS") {
    if (!("c_as" %in% names(par)))
      stop("Value for par$c_as missing!")
    if (!("alpha_as" %in% names(par)))
      stop("Value for par$alpha_as missing!")
  }
  
  # return parameters
  return(par)
}


.simulate_test = function(phy) {
  models = c("BM","OU","EB","JN","VG","NIG","AS","BMJN","BMVG","BMNIG","BMAS")
  par = list(sigma_bm=1,
             theta_ou=0.1,
             decay_eb=-0.1,
             lambda_jn=0.5,
             delta_jn=2,
             nu_vg=5,
             sigma_vg=2,
             alpha_nig=0.8,
             delta_nig=5,
             c_as=2,
             alpha_as=1.5,
             sigma_tip=0.1)
  
  phy = reorder(phy, "cladewise")
  x = list()
  for (m in models) {
    x[[m]] = rlevy(phy=phy, model=m, par=par)
    
    if (any(is.na(x[[m]]))) {
      stop(paste("simulations using",m,"generated NA values!"))
    }
  }
  
  return(x)
}

#####################
# solve for moments #
#####################

get_params_for_var = function(process_var=1,process_kurt=0,frac_of_var=0,halflife=0,decay=0,tip=0) {
  ret = list()
  ret$bm = list(sigma.bm=sqrt(process_var))
  ret$ou = get_params_OU(process_var,halflife)
  ret$eb = get_params_EB(process_var,decay)
  ret$jn = get_params_JN(process_var,process_kurt,1)
  ret$vg = get_params_VG(process_var,process_kurt,1)
  ret$nig = get_params_NIG(process_var,process_kurt,1)
  ret$bmjn = get_params_JN(process_var,process_kurt,frac_of_var)
  ret$bmvg = get_params_VG(process_var,process_kurt,frac_of_var)
  ret$bmnig = get_params_NIG(process_var,process_kurt,frac_of_var)
  ret$tip = tip
  return(ret)
}

get_params_JN = function(process_var,process_kurt,frac_of_var) {
  lambda = 3*frac_of_var^2/process_kurt
  delta2 = (process_kurt*process_var)/(3*frac_of_var)
  sigma2 = process_var-frac_of_var*process_var
  return(list(lambda.jn=lambda, delta.jn=sqrt(delta2), sigma.bm=sqrt(sigma2)))
}

get_params_VG = function(process_var,process_kurt,frac_of_var) {
  nu = process_kurt/(3*frac_of_var^2)
  tau2 = frac_of_var*process_var
  sigma2 = process_var-frac_of_var*process_var
  return(list(nu.vg=nu, sigma.vg=sqrt(tau2), sigma.bm=sqrt(sigma2)))
}

get_params_OU = function(process_var,halflife) {
  return(list(sigma.bm=sqrt(2*process_var*log(2)/halflife),alpha.ou=log(2)/halflife))
}

get_params_EB = function(process_var,decay) {
  # stationary variance of process is zero as t->Inf
  return(list(sigma.bm=sqrt(process_var),
              decay.eb=decay))
}

get_params_NIG = function(process_var,process_kurt,frac_of_var) {
  alpha = sqrt(3/(frac_of_var*process_kurt*process_var))
  delta = sqrt(3*frac_of_var*process_var/process_kurt)
  sigma2 = process_var - frac_of_var*process_var
  return(list(sigma.bm = sqrt(sigma2), alpha.nig = alpha, delta.nig = delta))
}

###################
# compute moments #
###################

get_moments_JN = function(lambda.jn, delta.jn, sigma.bm) {
  v = lambda.jn*delta.jn^2 + sigma.bm^2
  k = 3*lambda.jn*delta.jn^4/(sigma.bm^2+lambda.jn*delta.jn^2)^2
  return(list(var=v,kurt=k))
}

get_moments_VG = function(nu.vg, sigma.vg, sigma.bm) {
  v = sigma.bm^2+sigma.vg^2
  k = 3*nu.vg
  return(list(var = v, kurt = k))
}

get_moments_NIG = function(alpha.nig, delta.nig, sigma.bm) {
  v = sigma.bm^2 + delta.nig/alpha.nig
  k = 3*delta.nig/(alpha.nig*(sigma.bm^2*alpha.nig+delta.nig)^2)
  return(list(var = v, kurt = k))
}

#################
# Levy measures #
#################

JN_measure = function(x, lambda, delta, log = FALSE) {
  ret = log(lambda) + dnorm(x,sd=delta,log=TRUE)
  if (log) {
    return(ret)
  } else {
    return(exp(ret))
  }
}

VG_measure = function(x, nu, sigma, mu = 0 ,log=FALSE) {
  ret = -log(nu*abs(x)) + mu*x/sigma^2 - sqrt(2/nu-mu^2/sigma^2)/sigma*abs(x)
  if (log) {
    return(ret)
  } else {
    return(exp(ret))
  }
} 

NIG_measure = function(x, alpha, delta, beta = 0, log=FALSE) {
  ret = log(alpha) + log(delta) - log(pi) - log(abs(x)) + beta*x + log(besselK(alpha*abs(x),1))
  if (log) {
    return(ret)
  } else {
    return(exp(ret))
  }
}

#################
# rate of jumps #
#################

jump_rate = function(x, levy_measure, ...) {
  if (x == 0) {
    total = integrate(levy_measure, -Inf, Inf, ..., stop.on.error=FALSE)$value
  } else {
    
    negative = integrate(levy_measure, -Inf, -x, ..., stop.on.error=FALSE)$value
    positive = integrate(levy_measure, x, Inf, ..., stop.on.error=FALSE)$value
    total = negative + positive
  }
  return(total)
}

###########################
# moment from Levy measure#
###########################

moment_from_measure = function(k, levy_measure, ...) {
  integrate(function(y){y^k*levy_measure(y,...)},-Inf,Inf)
}


#########################
# brlen transformations #
#########################

OU.brlen = function(phy,theta=1e-6) {
  phy = reorder(phy,"postorder")
  n_tip = phy$Nnode + 1
  
  # get raw ages/times
  a = branching.times(phy)
  T = max(a)
  t = T - a
  
  # get OU-scaled ages/times
  t.ou = 1/(2*theta) * exp(-2*theta*(T-t)) * (1-exp(-2*theta*t))
  h.ou = 1/(2*theta) * exp(-2*theta*(T-T)) * (1-exp(-2*theta*T))
  a.ou = h.ou - t.ou
  a.ou = c( rep(0, n_tip), a.ou)
  
  # assign OU-scaled times to tree
  for (i in 1:nrow(phy$edge))
  {
    phy$edge.length[i] = a.ou[phy$edge[i,1]] - a.ou[phy$edge[i,2]]
  }
  return(phy)
}

EB.brlen = function(phy,r=1e-6) {
  phy = reorder(phy,"postorder")
  n_tip = length(phy$tip.label)
  
  # get raw ages/times
  a = branching.times(phy)
  T = max(a)
  t = T - a
  t = c(rep(T, n_tip), t)
  
  # assign EB-scaled times to tree
  for (i in 1:nrow(phy$edge))
  {
    t_pa = t[phy$edge[i,1]]
    t_ch = t[phy$edge[i,2]]
    dx = exp(r*t_ch) - exp(r*t_pa)
    phy$edge.length[i] = dx/r
  }
  
  return(phy)
}

####################
# parameter labels #
####################

# parameter labeling
format_params = function(p, m) {
  param_names = c("sigma_bm",
                  "lambda_jn",
                  "delta_jn",
                  "sigma_vg",
                  "nu_vg",
                  "alpha_nig",
                  "delta_nig",
                  "alpha_ou",
                  "decay_eb",
                  "sigma_tip")
  n_param = length(param_names)
  x = rep(0, n_param)
  names(x) = param_names
  
  x[n_param] = p[length(p)]
  if (m=="BM") {
    x[1]=p[1]
  } else if (m=="BMJN") {
    x[1]=p[1]
    x[2]=p[2]
    x[3]=p[3]
  } else if (m=="BMVG") {
    x[1]=p[1]
    x[4]=p[2]
    x[5]=p[3]
  } else if (m=="BMNIG") {
    x[1]=p[1]
    x[6]=p[2]
    x[7]=p[3]
  } else if (m=="JN") {
    x[2]=p[1]
    x[3]=p[2]
  } else if (m=="VG") {
    x[4]=p[1]
    x[5]=p[2]
  } else if (m=="NIG") {
    x[6]=p[1]
    x[7]=p[2]
  } else if (m=="OU") {
    x[1]=p[1]
    x[8]=p[2]
  } else if (m=="EB") {
    x[1]=p[1]
    x[9]=p[2]
  }
  return(x)
}


#################
# data cleaning #
#################

drop.outlier = function(phy,dat,n=1,drop_zero=T,verbose=F)
{
  
  # expects dat to be a named vector
  if (!is.vector(dat))
  {
    stop("ERROR: dat is not a vector")
  }
  if (is.null(names(dat)))
  {
    stop("ERROR: dat does not have names")
  }
  
  to.drop = c()
  to.drop.zero = c()
  to.drop.outlier = c()
  
  if (drop_zero)
  {
    td = drop.zero(phy,dat)
    phy = td$phy
    dat = td$dat
    if (length(td$dropped) > 0)
    {
      to.drop.zero = td$dropped
    }
  }
  
  phy = reorder(phy,'postorder')
  dat = dat[phy$tip.label]
  
  # only drop outliers if n > 0 
  if (n > 0) {
    contrast = pic(dat,phy)
    nodes = order(abs(contrast),decreasing=TRUE)[1:n]
    nodes = length(phy$tip.label)+nodes
    for (node in nodes)
    {
      clade = extract.clade(phy,node)
      to.drop.outlier = clade$tip.label
    }
  }
  to.drop = unique(c(to.drop.outlier, to.drop.zero))
  
  #to.drop = unique(to.drop)
  #cat("to.drop\n")
  #print(to.drop)
  #cat("phy$tip.labels\n")
  #print(phy$tip.label)
  #print(c(length(to.drop), length(phy$tip.label)))
  
  if (length(to.drop) != length(phy$tip.label))
  {
    newPhy = drop.tip(phy, to.drop)
    newDat = dat[newPhy$tip.label]
  } else {
    newPhy = phy
    newDat = dat[newPhy$tip.label]
  }
  if (verbose) {
    cat("Dropped taxa\n")
    print(to.drop)
    cat("Dropped contrasts\n")
    print(contrast[order(abs(contrast),decreasing=TRUE)[1:n]])
    cat("Dropped taxa (drop_zero)    =",length(to.drop.zero),"\n")
    cat("Dropped taxa (drop_outlier) =",length(to.drop.outlier),"\n")
    cat("\n")
  }
  return(list(phy=newPhy,dat=newDat))
}

drop.zero = function(phy,dat,eps=0) {
  #expects dat to be a named vector
  if (!is.vector(dat)) {
    stop("ERROR: dat is not a vector")
  }
  if (is.null(names(dat))) {
    stop("ERROR: dat does not have names")
  } 
  phy = reorder(phy,"postorder")
  dat = dat[phy$tip.label]
  bad_nodes = c()
  pseudo_obs = rep(0,Nnode(phy,internal.only=FALSE))
  for (i in 1:length(phy$tip.label)) {
    pseudo_obs[i] = dat[i]
  }
  for (i in seq(1,nrow(phy$edge),2)) {
    cur_node = phy$edge[i,1]
    cur_left = phy$edge[i,2]
    cur_right = phy$edge[i+1,2]
    t_l = phy$edge.length[i]
    t_r = phy$edge.length[i+1]
    contrast = pseudo_obs[cur_left]-pseudo_obs[cur_right]
    if (!is.nan(contrast) && abs(contrast) <= eps) {
      bad_nodes = c(bad_nodes,cur_node)
    }
    pseudo_obs[cur_node] = t_r*pseudo_obs[cur_left]+t_l*pseudo_obs[cur_right]
    pseudo_obs[cur_node] = pseudo_obs[cur_node]/(t_l+t_r)
    
  }
  to.drop = c()
  for (node in bad_nodes) {
    clade = extract.clade(phy,node)
    to.drop = c(to.drop,clade$tip.label)
  }
  to.drop = unique(to.drop)
  newPhy = drop.tip(phy,to.drop)
  newDat = dat[newPhy$tip.label]
  return(list(phy=newPhy,dat=newDat,bad_nodes=bad_nodes,dropped=to.drop))
}


### plot PIC outliers for BM component of BM+LP

pic_outliers = function(phy,dat,sigma_bm) {
  cur_pic = pic(dat,phy)
  z_vals = cur_pic
  p_vals = pnorm(abs(z_vals),sd=sigma_bm,lower.tail=FALSE)
  return(-log10(p_vals))
}

get_p_cols = function(p_vals,phy) {
  normed = p_vals/max(p_vals)
  cols_rgb = colorRamp(c("white","red"),bias=100)(normed)
  cols_hex = apply(cols_rgb,1,function(x){rgb(x[1],x[2],x[3],maxColorValue=255)})
  names(cols_hex) = names(p_vals)
  return(cols_hex)
}

plot_jumps = function(x,cutoff=-log10(0.05),cex=.5,adj=.5,main="") {
  
  if (is.null(x$phy)) stop("x does not contain phy object!")
  if (is.null(x$dat)) stop("x does not contain dat object!")
  if (is.null(x$params)) stop("x does not contain params object!")
  if (!("sigma_bm" %in% names(x$params))) stop("x does not contain sigma_bm parameter!")
  
  sigma = x$params["sigma_bm"]
  dat_p_vals = pic_outliers(x$phy, x$dat, sigma)
  max_p_val = max(dat_p_vals[dat_p_vals != Inf])
  if (any(dat_p_vals==Inf)) {
    inf_idx = dat_p_vals==Inf
    dat_p_vals[inf_idx] = 1.5 * max_p_val
  }
  filt_p_vals = dat_p_vals[dat_p_vals>=cutoff]
  
  p_cols = get_p_cols(filt_p_vals)
  plot(x$phy,adj=adj,cex=cex,main=main)
  nodelabels(pch=16,node = as.numeric(names(p_cols)), col=p_cols,frame="circle")
  invisible(dat_p_vals)
}