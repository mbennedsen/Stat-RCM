# Stat-RCM
Implementation of Stat-RCM from the paper Bennedsen, Hillebrand, and Koopman (2025), “A Statistical Reduced Complexity Climate Model for Probabilistic Analyses and Projections”, Journal of Climate, Volume 38, Issue 21, pp. 6329-6350.

This code and data will reproduce the main output from the paper. The code may be used, distributed, and changed freely.

The folder ‘Files’ contain pre-estimated parameters and bootstrap samples (for use in calculating standard errors and assessing uncertainty in the estimated ECS) for Stat-RCM with the sqrt+log forcings equation (found to provide the lowest BIC, see the main paper) applied to the full data period 1959-2022.

The ‘main.m’ script may be run using different data (e.g. updated data over a longer sample period), but the estimations, especially the bootstrapping, can be time-consuming. If bootstrap standard errors are required, it is recommended to parallelize the bootstrap script “bootstrap_stderrs_sqrtlog.m” using a parfor-loop.

The ‘projections.m’ script may be run to produce projections using the estimated Stat-RCM model and a future SSP scenario. Other scenarios may be considered by altering the code slightly.
