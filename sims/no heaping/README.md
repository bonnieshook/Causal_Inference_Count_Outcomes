Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women

R code for simulations without data heaping
These programs implement the simulation study described in Shook-Sa et al, "Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women" in R

The simulations are implemented in the following order:

1. Generate and analyze data for the simulation study for each sample size and distribution of interest using, for example, the 02_Poisson.R program, called with the 02_Run_Sims_Poisson_800.sh shell script. This process was repeated for each of the other three distributions (negative binomial, ZIP, and ZINB).
   Note that these programs call the 00_Estimators_11.03.21.R program, which contains functions for computing each estimator of interest along with corresponding 95% confidence intervals.
2. Combine the results from all iterations of the simulation using the 06_compile_sims.R program, called with the shell script 06_Compile_Sims_all.sh. This program creates a final dataset with one row for each iteration of the simulations, for each distribution of interest.
3. Calculate summary measures for each distribution using 07_combine_sim_tables.R.

Developed by: Bonnie Shook-Sa
