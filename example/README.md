Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women

R code for analyzing a single dataset without or with data heaping ('noheap.csv' or 'heaping.csv')


### Without data heaping
The dataset 'noheap.csv' contains the following variables:

part.id = participant id (1-n),
agevis10 = continuous covariate,
druguse = binary covariate,
sxex = binary covariate,
inc = binary exposure,
SP = count outcome (no heaping)

The program 01_Analysis_no_heaping.R demonstrates how to analyze the example dataset using the IPTW, parametric g-formula, and doubly robust estimators demonstrated in the paper, assuming that there is no data heaping.
This program calls the 00_Estimators_no_Heaping_10.23.23.R program, which contains the functions for each estimator. Both programs are annotated with the locations that should be updated if additional covariates are added or removed
from the weight or outcome models.


### With data heaping
The dataset 'heaping.csv' contains the following variables:

part.id = participant id (1-n),
income = continuous covariate,
inc = binary exposure,
CIG.heap = count outcome (with heaping)

The program 02_Analysis_heaping.R demonstrates how to analyze the example dataset using the IPTW, parametric g-formula, and doubly robust estimators demonstrated in the paper, assuming that there is data heaping. All four estimators (3 HCAR and 1 IH) are demonstrated.
This program calls the 00_Estimators_Heaping_10.23.23.R program, which contains the functions for each estimator. Both programs are annotated with the locations that should be updated if additional covariates are added or removed
from the weight or outcome models.

Developed by: Bonnie Shook-Sa
