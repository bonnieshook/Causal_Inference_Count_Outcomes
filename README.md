# Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women

### Bonnie E. Shook-Sa, Michael G. Hudgens, Andrea K. Knittel, Andrew Edmonds, Catalina Ramirez, Stephen R. Cole, Mardge Cohen, Adebola Adedimeji, Tonya Taylor, Katherine G. Michel, Andrea Kovacs, Jennifer Cohen, Jessica Donohue, Antonina Foster, Margaret A. Fischl, Dustin Long, Adaora A. Adimora

**Citation**: Shook-Sa BE, Hudgens MG, Knittel AK, Edmonds A, Ramirez C,  Cole SR, Cohen M, Adedimeji A, Taylor T, Michel KG, Kovacs A, Cohen J, Donohue J, Foster A, Fischl MA, Long D, Adimora AA. "Exposure Effects on Count Outcomes with Observational Data, with Application to Incarcerated Women." *arXiv* arXiv:2202.01650.
[![arXiv](https://img.shields.io/badge/arXiv-2202.01650-b31b1b.svg)](https://arxiv.org/abs/2202.01650)
--------------------------------

## Abstract

Causal inference methods can be applied to estimate the effect of a point exposure or treatment on an outcome of interest using data from observational studies. For example, in the Women's Interagency HIV Study, it is of interest to understand the effects of incarceration on the number of sexual partners and the number of cigarettes smoked after incarceration. In settings like this where the outcome is a count, the estimand is often the causal mean ratio, i.e., the ratio of the counterfactual mean count under exposure to the counterfactual mean count under no exposure. This paper considers estimators of the causal mean ratio based on inverse probability of treatment weights, the parametric g-formula, and doubly robust estimation, each of which can account for overdispersion, zero-inflation, and heaping in the measured outcome. Methods are compared in simulations and are applied to data from the Women's Interagency HIV Study.

--------------------------------

## File Manifesto

### example
The `example/` path contains example datasets with and without data heaping and code for analyzing the data using the methods presented in the paper.

### sims
The `sims/` path contains R programs to replicate the simulation study presented in the paper.
