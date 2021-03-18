# echinoid_body_size_2021

## Introduction
Repository for data files and code associated with a macroevolutionary analysis of body size in echinoids (sea urchin, heart urchin and sand dollars). This study is published in: Mongiardino Koch (2021). Exploring adaptive landscapes across deep time: a case study using echinoid body size. *Evolution*, in press.

## Description
The main file in this repository is `bodysize_analyses.R` which can be used to replicate all stages of model fitting, output some statistics as well as some summary plots. This includes the comparison of a wide range of uniform macroevolutionary models against multi-peak Ornstein-Uhlenbeck processes that capture evolutionary dynamics across adaptive landscapes. The code also makes direct us of the [extendedSurface](https://github.com/mongiardino/extendedSurface) algorithm (developed for this study), which was found to outperform other approaches of constraining the number and location of adaptive regimes shifts without the need for these to be specified *a priori*.

Other files deposited here include:

- `echinoid_tree.tre`: the maximum clade credibility (mcc) topology from the total-evidence dated analysis of Mongiardino Koch & Thompson (2020). A Total-Evidence Dated Phylogeny of Echinoidea Combining Phylogenomic and Paleontological Data, *Systematic Biology*, syaa069, https://doi.org/10.1093/sysbio/syaa069.

- `random_posterior_trees.tre`: a random selection of 20 topologies from the same study.

- `bodysize_dataset.txt`: the dataset with body sizes and error values for tips in the analysis.

- `accesory_funtions.R`: some functions needed for the analyses.
