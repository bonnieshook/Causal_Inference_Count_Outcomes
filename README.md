# Causal_Inference_Count_Outcomes

This repository contains simulation code to replicate the findings in "Estimating Point Exposure Effects on Count Outcomes with Observational Data" by Bonnie E. Shook-Sa, Michael G. Hudgens, et al.

*00_Estimators_11.03.21 contains functions for point and variance estimation in the absence of data heaping
*00_Estimators_Heaping_01.09.22 contains functions for point and variance estimation that account for data heaping
*02_Poisson_11.03.21, 03_NBin_11.03.21, 04_ZIP_11.03.21, and 05_ZINB_11.03.21 conduct the simulations for each of the four parametric outcomes considered in the manuscript in the absence of data heaping
*06_combine_sim_tables_11.03.21 combines the results of the simulations in the absence of data heaping
*10_Heaping_corr_sp_01.09.22 and 10_Heaping_inc_sp_01.09.22 conduct the simulations with data heaping present under correct and incorrect model specification, respectively
*11_combine_sim_tables_heap combines the results of the simulations in the presence of data heaping
