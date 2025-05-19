# Sampling_Conditional_Vines

This is the code corresponding to the paper "Sampling from Conditional Distributions of Simplified Vines"
by Ariane Hanebeck, Özge Sahin, Petra Havlickova, Claudia Czado
published in Statistics and Computing.

We provide the code for two examples - a univariate and a bivariate conditional distribution.

The code for sampling from a conditional distribution of a specified vine is given in:
- Example1dim.R for the univariate case,
- Example2dim.R for the bivariate case.

These are the two files that the user needs to work with.

The other files contain the internal functions:
- Sample_Function.R contains the function, which starts sampling with Stan,
- Transformation.R allows for a transformation from vines defined in rvinecopulib to vines defined in VineCopula.

There are two Stan-files:
- STAN.stan contains the Stan-code for univariate sampling,
- STAN2.stan contains the Stan-code for multivariate sampling.