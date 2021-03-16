########################################################################################################################

#This script contains code to perform macroevolutionary analyses of echhinoid body size.
#Datasets required to run these analyses are also provided.

#Code by Nicolás Mongiardino Koch during 2020 covid-19 lockdown, to
#perform research in Mongiardino Koch N - Exploring adaptive landscapes across
#deep time: a case study using echinoid body size

########################################################################################################################

##Setup-------------------------------------------------------------------------------------------------
rm(list=ls())

#Change working directory to location of files that are provided with this R
#code file
setwd('')

#install required packages
packages <- c('ape','surface','OUwie','geiger','BBMV','phytools','tidyverse', 'devtools')
new_packages <- packages[!packages %in% installed.packages()[,'Package']]
if(length(new_packages)) { install.packages(new_packages) }

#install the package developed for this analysis: extendedSurface
devtools::install_github("mongiardino/extendedSurface")

library(ape)
library(surface)
library(OUwie)
library(extendedSurface)
library(geiger)
library(BBMV)
library(phytools)
library(tidyverse)

#load accompanying file including a few useful functions as well as the entire
#pulseR package which fails to install directly from github
source('accesory_funtions.R')


##Analyses-------------------------------------------------------------------------------------------------
#load data
size_data = read.csv('bodysize_dataset.txt')
number_of_specimens = size_data[,3]
size = size_data[,4]
error = size_data[,5]
names(number_of_specimens) = names(size) = names(error) = size_data[,2]

#load trees
mcc = read.tree('echinoid_tree.tre')
trees_posterior = read.tree('random_posterior_trees.tre')
for(i in 1:length(trees_posterior)) {
  trees_posterior[[i]] = drop.tip(trees_posterior[[i]], which(trees_posterior[[i]]$tip.label %not in% mcc$tip.label))
}

#remove all fossils
root_tip = round(dist.nodes(mcc)[(length(mcc$tip.label)+1), (1:length(mcc$tip.label))], 1)
timespan = round(max(nodeHeights(mcc)), 1)
fossil = mcc$tip.label[as.numeric(which(root_tip < timespan))]
tree_extant_A = drop.tip(mcc, fossil)
tree_extant_A = force.ultrametric(tree_extant_A, method = 'extend')

#remove fossil lineages 
other_extant_lineages = c('Arachnoididae', 'Carnarechinus_clypeatus', 'Echinoneus_cyclostomus', 'Glyptocidaris_crenularis', 
                          'Goniocidaris_tubaria', 'Kamptosoma_asterias', 'Maretia_planulata', 'Paraphormosoma_alternans', 
                          'Spatangus_purpureus', 'Coelopleurus', 'Corystus', 'Hemiaster_bufo', 'Macropneustes', 'Pericosmus', 
                          'Salenia_petalifera', 'Temnopleurus_toreumaticus', 'Prenaster_alpinus', 'Scutella_subrotunda')

tree_extant_B = drop.tip(mcc, fossil[which(fossil %not in% other_extant_lineages)])
tree_extant_B = force.ultrametric(tree_extant_B, method = 'extend')

#Uncomment one of the next lines to choose which tree should be used, i.e. the
#mcc, a sample of 20 posterior topologies, or two types of extant trees,
#corresponding to the extant terminals (tree_extant_A) or those plus the sampled
#fossils from extant clades stretched to the present (tree_extant_B)
tree = mcc
#tree = trees_posterior
#tree = tree_extant_A
#tree = tree_extant_B

#subsample size data to those present in the topology chosen
size_to_use = size[which(names(size) %in% tree$tip.label)]
error_to_use = error[which(names(error) %in% tree$tip.label)]

#Prepare datasets to hold AICc values for all models
models = c('BM', 'OU', 'EB', 'Trend', 'White', 'BBM', 'FPK', 'BBMV', 'Pulsed', 'OUM', 'BMS', 'OUMVA', 'OUMVAZ')
AICC = matrix(ncol = length(models), nrow = ifelse(class(tree) == 'multiPhylo', length(tree), 1))
colnames(AICC) = models

#fit all 13 models of macroevolution to either one tree or all trees in the
#multiphylo object chosen
for(i in 1:length(tree)) {
  if(class(tree) == 'multiPhylo') {
    this_tree = tree[[i]]
  } else {
    this_tree = tree
  }
  #BM
  BM = fitContinuous(this_tree, size_to_use, SE = error_to_use, model = 'BM')
  AICC[i,1] = aic_to_aicc(aic = (-2*BM$opt$lnL) + (2*(BM$opt$k-1)), param = BM$opt$k-1, tree = this_tree)
  
  #OU
  OU = fitContinuous(this_tree, size_to_use, SE = error_to_use, model = 'OU')
  AICC[i,2] = aic_to_aicc(aic = (-2*OU$opt$lnL) + (2*(OU$opt$k-1)), param = OU$opt$k-1, tree = this_tree)
  
  #EB
  EB = fitContinuous(this_tree, size_to_use, SE = error_to_use, model = 'EB')
  AICC[i,3] = aic_to_aicc(aic = (-2*EB$opt$lnL) + (2*(EB$opt$k-1)), param = EB$opt$k-1, tree = this_tree)
  
  #Trend
  trend = fitContinuous(this_tree, size_to_use, SE = error_to_use, model = 'trend')
  AICC[i,4] = aic_to_aicc(aic = (-2*trend$opt$lnL) + (2*(trend$opt$k-1)), param = trend$opt$k-1, tree = this_tree)
  
  #White noise
  whitenoise = fitContinuous(this_tree, size_to_use, SE = error_to_use, model = 'white')
  AICC[i,5] = aic_to_aicc(aic = (-2*whitenoise$opt$lnL) + (2*(whitenoise$opt$k-1)), param = whitenoise$opt$k-1, tree = this_tree)
  
  #BBM, FPK and BBMV do not incorporate error estimates directly, rather they
  #use vectors of plausible observations given error values. The following
  #genertates these for all taxa
  size_list = list()
  for(j in 1:length(size_to_use)) {
    size_list[[j]] = rnorm(n = 10, mean = size_to_use[j], 
                           sd = sqrt(error_to_use[j])*sqrt(number_of_specimens[j]))
  }
  names(size_list) = names(size_to_use)
  bounds = c(min(unlist(size_list)), max(unlist(size_list)))
  
  #Bounded Brownian motion
  ll_BBMV0 = lnL_BBMV(this_tree, size_list, Npts = 100, bounds = bounds, a = 0, b = 0, c = 0)
  fitBBM = find.mle_FPK(model = ll_BBMV0, method = 'L-BFGS-B', safe = T)
  AICC[i,6] = aic_to_aicc(aic = ((-2*fitBBM$lnL) + ((fitBBM$k+1)*2)), param = (fitBBM$k+1), tree = this_tree)
  
  #Fokker-Plank-Kolmogorov
  ll_FPK4 = lnL_FPK(this_tree, size_list, Npts = 100, a = NULL, b = NULL, c = NULL)
  fitFPK4 = find.mle_FPK(model = ll_FPK4)
  AICC[i,7] = aic_to_aicc(aic = ((-2*fitFPK4$lnL) + (fitFPK4$k*2)), param = fitFPK4$k, tree = this_tree)
  
  #BBMV
  ll_BBMV4 = lnL_BBMV(this_tree, size_to_use, Npts = 100, bounds = bounds, a = NULL, b = NULL, c = NULL)
  ll_BBMV2 = lnL_BBMV(this_tree, size_to_use, Npts = 100, bounds = bounds, a = 0, b = NULL, c = NULL)
  ll_BBMV1 = lnL_BBMV(this_tree, size_to_use, Npts = 100, bounds = bounds, a = 0, b = 0, c = NULL)
  fit4b = find.mle_FPK(model = ll_BBMV4)
  fit2b = find.mle_FPK(model = ll_BBMV2)
  fit1b = find.mle_FPK(model = ll_BBMV1)
  BBMV = c(aic_to_aicc(aic = ((-2*fit4b$lnL) + (fit4b$k*2)), param = fit4b$k, tree = this_tree), 
           aic_to_aicc(aic = ((-2*fit2b$lnL) + (fit2b$k*2)), param = fit2b$k, tree = this_tree), 
           aic_to_aicc(aic = ((-2*fit1b$lnL) + (fit1b$k*2)), param = fit1b$k, tree = this_tree))
  names(BBMV) = c('3param', '2param', '1param')
  AICC[i,8] = as.numeric(BBMV[which(BBMV == min(BBMV))])
  cat('From BBMV the winner is', names(BBMV[which(BBMV == min(BBMV))]), '\n')
  
  #Pulsed
  fitJN = fit_reml_levy(this_tree, size_to_use, model = 'JN', sigma_tip = T)
  fitNIG = fit_reml_levy(this_tree, size_to_use, model = 'NIG', sigma_tip = T)
  fitBMJN = fit_reml_levy(this_tree, size_to_use, model = 'BMJN', sigma_tip = T)
  fitBMNIG = fit_reml_levy(this_tree, size_to_use, model = 'BMNIG', sigma_tip = T)
  pulsed = c(aic_to_aicc(fitJN$AIC, fitJN$n_params, this_tree), 
             aic_to_aicc(fitNIG$AIC, fitNIG$n_params, this_tree), 
             aic_to_aicc(fitBMJN$AIC, fitBMJN$n_params, this_tree), 
             aic_to_aicc(fitBMNIG$AIC, fitBMNIG$n_params, this_tree))
  names(pulsed) = c('JN', 'NIG', 'BMJN', 'BMNIG')
  AICC[i,9] = as.numeric(pulsed[which(pulsed == min(pulsed))])
  cat('From pulsed the winner is', names(pulsed[which(pulsed == min(pulsed))]), '\n')
  
  #Prepare data for SURFACE
  size_surface = data.frame(size_to_use)
  tree_surface = nameNodes(this_tree)
  data_surface = convertTreeData(tree_surface, size_surface)
  
  #SURFACE
  fwd.surface = surfaceForward(data_surface[[1]], data_surface[[2]], verbose = T)
  bwd.surface = surfaceBackward(data_surface[[1]], data_surface[[2]], starting_model = fwd.surface[[length(fwd.surface)]], verbose = T, aic_threshold = 0)
  AICC[i,10] = aic_to_aicc(bwd.surface[[length(bwd.surface)]]$aic, bwd.surface[[length(bwd.surface)]]$n_regimes[2]+2, this_tree)
  
  #Prepare data for OUwie
  size_ouwie = data.frame(rownames(size_surface), NA, size_surface$size_to_use, error_to_use)
  for(j in 1:nrow(size_ouwie)) {
    size_ouwie[j,2] = as.character(unlist(attributes(bwd.surface[[length(bwd.surface)]]$fit$size_to_use)[6])[which(as.character(unlist(attributes(bwd.surface[[length(bwd.surface)]]$fit$size_to_use)[15])) == as.character(size_ouwie[j,1]))])
  }
  tree_ouwie = tree_surface
  for(j in 1:length(tree_surface$node.label)) {
    tree_ouwie$node.label[j] = as.character(unlist(attributes(bwd.surface[[length(bwd.surface)]]$fit$size_to_use)[6]))[which(as.character(unlist(attributes(bwd.surface[[length(bwd.surface)]]$fit$size_to_use)[15])) == tree_surface$node.label[j])]
  }
  
  #multi Brownian motion
  BMS = OUwie(tree_ouwie, size_ouwie, model = 'BMS', simmap.tree = F, root.age = tree_ouwie$root.time, scaleHeight = F, 
              get.root.theta = F, clade = NULL, mserr = 'known', starting.vals = NULL)
  AICC[i,11] = aic_to_aicc((-2*(BMS$loglik)) + (2*(BMS$param.count-1)), BMS$param.count-1, this_tree)
  
  #OUMVA and OUMVAZ (depending on whether trees are ultrametric)
  if(is.ultrametric(this_tree)) {
    extended = surfaceExtended(bwd.surface, size_surface, tree_surface, error = error_to_use, models = "OUMVA", 
                               limit = 2, plot = T, fwd_surface = fwd.surface)
  } else {
    extended = surfaceExtended(bwd.surface, size_surface, tree_surface, error = error_to_use, models = c("OUMVA", "OUMVAZ"), 
                               limit = 2, plot = T, fwd_surface = fwd.surface)
    AICC[i,13] = min(extended$summary[endsWith(as.character(extended$summary[,1]),'OUMVAZ'),3], na.rm = T)
  }
  
  AICC[i,12] = min(extended$summary[endsWith(as.character(extended$summary[,1]),'OUMVA'),3], na.rm = T)
  
  #if the tree is a multiphylo object, save the OUM/OUMVA/OUMVAZ models for each
  #tree to plot some results
  if(class(tree) == 'multiPhylo') {
    assign(paste0('OUM_', i), bwd.surface)
    assign(paste0('OUMVA_', i), extended$best_OUMVA)
    assign(paste0('OUMVAZ_', i), extended$best_OUMVAZ)
  } else {
    break
  }
}

##Results---------------------------------------------------------------------------------------------------
if(!class(tree) == 'multiPhylo') {
  #rank models based on fit
  aicw(AICC[,which(!is.na(AICC))])
} else {
  aic = tibble(aic = c(AICC), 
               model = rep(models, each = length(tree)), 
               tree = rep(1:length(tree), length(models))) %>% arrange(tree)
  
  weights = c()
  for(i in 1:length(tree)) {
    treex = unlist(aic[which(aic$tree == i),1], use.names = F)
    names(treex) = unlist(aic[which(aic$tree == i),2], use.names = F)
    weights = c(weights, aicw(treex)[]$aicweights)
  }
  aic = cbind(aic, weights)
  
  medians = vector(length = length(models))
  names(medians) = unique(aic$model)
  for(i in 1:length(models)) {
    medians[i] = median(unlist(aic[which(aic$model == unique(aic$model)[i]),1], use.names = F))
  }
  
  aic$model = factor(aic$model, levels = names(sort(medians, decreasing = T)))
  medians = data.frame(aic = sort(medians, decreasing = T), 
                       model = factor(levels(aic$model), levels = levels(aic$model)))
  
  #Plot model fitting results across trees(Fig. 2B of manuscript)
  ggplot(aic, aes(x = model, y = aic)) + geom_line(aes(group = tree), color = '#bebebe', size = 0.5) + 
    geom_point(aes(fill = weights), size = 5, shape = 21, color = 'black', alpha = 0.6) + 
    geom_line(data=medians, aes(x = model, y = aic, group = NA), size = 1, color = '#231f20') + 
    scale_fill_gradient(high = '#de2d26', low = 'white') + theme(legend.position = "none")
  
  #Inset of Fig.2B with number of regimes across different types of multi-OU
  #parameterizations
  to.plot = data.frame(optima = NA, model = rep(c('OUM', 'OUMVA', 'OUMVAZ'), each = 20))
  for(i in 1:20) {
    a = get(paste0('OUM_', i))
    to.plot[i,1] = as.numeric(a[[length(a)]]$n_regimes[2])
    a = get(paste0('OUMVA_', i))
    to.plot[(i+20),1] = ncol(a$solution)
    a = get(paste0('OUMVAZ_', i))
    to.plot[(i+40),1] = ncol(a$solution)
  }
  
  X = min(to.plot[,1]):max(to.plot[,1])
  Y = unique(to.plot[,2])
  data = expand.grid(X = X, Y = Y)
  Z = rep(NA, nrow(data))
  data = cbind(data, Z)
  
  for(i in 1:nrow(data)) {
    a = to.plot[which(to.plot[,2] == data$Y[i]),]
    data$Z[i] = length(which(a$optima == data$X[i]))/20
  }
  
  ggplot(data, aes(X, Y, fill= Z)) + geom_tile() + theme_classic() + 
    scale_fill_gradient(low = 'white', high = '#ff8201') + coord_flip()
}