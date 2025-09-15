# Additive genetic variance for fitness in _Brassica rapa_

### Accepted: The capacity for adaptation to climate warming in a naturalized annual plant (_Brassica rapa_)
### Authors: Cameron So, Karl Grieshop, Sydney Rotman, Arthur Weis

This project has been accepted in the journal Evolution. 
A pre-print was released in 2021 on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.01.510426v1). 
The data and code for this project will also be soon archived on [Dryad](datadryad.org). 

### Abstract

The persistence of a declining population in the face of environmental change may depend on how fast natural selection restores fitness, a process called “evolutionary rescue”. In turn, evolutionary rescue depends on a population’s adaptive potential. Fisher’s theorem states that a population’s adaptive potential equals the additive genetic variance for fitness (_V<sub>A</sub>_(_W_)) divided by mean fitness Embedded Image. Both the numerator and denominator of this rate can differ across environments even when holding allele frequencies constant. However, little is known about how these rates change in wild populations during adaptation, including changes in additive and dominance variance. We assessed the change in adaptive potential and dominance variance in fitness (_V<sub>D</sub>_(_W_)) for a Québec population of wild mustard (_Brassica rapa_) under climate warming. We also assessed adaptive constraints that could arise from negative genetic correlations across environments. We grew a pedigreed population of 7000 plants under ambient and heated (+4°C) temperatures and estimated the change in Embedded Image, _V<sub>A</sub>_(_W_), _V<sub>D</sub>_(_W_), and the cross-environment genetic correlations (_r<sub>A</sub>_). As predicted, estimates of _V<sub>A</sub>_(_W_) and adaptive potentials were higher under heated conditions but non-significantly so. This is perhaps because, surprisingly, plants exposed to a warmer climate exhibited greater Embedded Image. Nevertheless, increased fitness in the warmer environment suggests a plasticity-based short-term potential for adaptation, and that weak but non-significant genetic correlations across environments will enable slow on-going adaptation to warming. Overall, this population of _B. rapa_ harbours existing genetic architecture to persist under warmer temperatures through pre-adaptation but not through evolutionary rescue

It includes quantitative genetic models to estimate additive genetic variance for fitness components and growth traits, using MCMCglmm. 

### How to use this code

R is the main statistical software.

Main scripts:
- 01_Data_exploration.R - Data processing and exploratory graphs
- 02_MCMC_data_prep.R - Full data processing step, with overlaps from 01_Data_exploration.R
- 03_MCMC_ambient.R - Estimates quantitative genetic parameters (Va, Vd, Vm, h2) and trait means for plants in the ambient treatment
- 04_MCMC_heated.R - Estimates quantitative genetic parameters (Va, Vd, Vm, h2) and trait meansfor plants in the heated treatment
- 05_MCMC_genetic_correlations.R - Estimates additive genetic variance for plasticity and additive genetic correlations
- 06_MCMC_Climate_analysis.R - Extracts and plots climatic data

Supplementary scripts:
- S1_Frequentist_trait_mean_analysis.R - Estimates trait means and treatment effects under Frequentist framework
- S2_Flowering_time_analysis.R - Estimates flowering times and treatment effects
- S3_MCMC_full_models.R - Attempts to estimate quantitative genetic parameters in full model, including ambient and heated data
- S4_lme4_Vg_broad.R - Estimates broad sense genetic variance under Frequentist framework
- S5_MCMC_Vg_broad.R - Estimates broad sense genetic variance under Bayesisan framework

### Acknowledgements
 
* **Frank Shaw** - *R Code* - [Email](mailto:fshaw314@gmail.com?subject=[GitHub]%20Source%20Han%20Sans)
* **Charlie Geyer** - *R Code* - [Email](mailto:geyer@umn.edu?subject=[GitHub]%20Source%20Han%20Sans)
* **Pierre de Villemereuil** - *MCMCglmm assistance* - [Email](mailto:pierre.de-villemereuil@mnhn.fr?subject=[GitHub]%20Source%20Han%20Sans)
* **Matthew Wolak** - *nadiv assistance* - [Email](mailto:terps@auburn.edu?subject=[GitHub]%20Source%20Han%20Sans)
