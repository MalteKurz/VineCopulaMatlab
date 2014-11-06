# PairCopulaParameterBounds

Getting the parameter bounds for pair-copulas

---
## Purpose
The function gives back bounds for the parameter values. Those can
        be set globally for the VineCPP toolbox by using the GUI
        SetGlobalParameterBounds. The bounds are used for checking for
        correctly specified pair-copulas and also as bounds for the
        estimation of the parameters.


## Usage
               [lb,ub] = PairCopulaParameterBounds(families)


## Inputs
       families        = Either the copula family or a cell array or
                         numerically coded vector of copula families.


## Outputs
       lb              = If the input variable families is just one single
                         copula, then it is the vector of lower bounds of
                         the corresponding copula parameter(s). If the input
                         variable families is a cell array of copula
                         families, then it is a vector of lower bounds for
                         all the copulas being specified in the families
                         cell array.
       ub              = If the input variable families is just one single
                         copula, then it is the vector of upper bounds of
                         the corresponding copula parameter(s). If the
                         input variable families is a cell array of copula
                         families, then it is a vector of upper bounds for
                         all the copulas being specified in the families
                         cell array.
