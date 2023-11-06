Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women

R code for simulations with data heaping
These programs implement the simulation study described in Shook-Sa et al, "Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women" in R

The simulations are implemented in the following order:

1. Generate and analyze data for the heaping simulation study for each sample size, scenario (IH and HCAR), and model specification (correct or incorrect) of interest using, for example, the 01_Heaping_corr_sp_IH.R program, called with the 01_Run_Sims_Heap_Corr_Sp_IH_800.sh shell script. This process was repeated with the 02, 04, and 05 programs for correctly specified HCAR and incorrectly specified IH and HCAR scenarios, respectively.
   Note that these programs call the 00_Estimators_11.03.21.R and 00_Estimators_Heaping_07.17.23.R programs, which contains functions for computing each naive (nonheaping) and heaping estimator of interest along with corresponding 95% confidence intervals.
2. Combine the results from all iterations of the simulation using the 10_compile_sims_heap_exclude.R program, called with the shell script 10_Compile_Sims_heap_800.sh. This program creates a final dataset with one row for each iteration of the simulations, for each scenario of interest.
3. Calculate summary measures for each scenario using 11_combine_sim_tables_heap.R.

Developed by: Bonnie Shook-Sa
