# PairCopulaKendallTau

Transforming the parameter of a pair-copula into the value of Kendall's tau

---
## Purpose
The function computes the value of Kendall's tau given the
        parameter of the specific two-dimensional copula. Possible copula
        families:
        
           1   AMH
           7   Clayton
           8   FGM
           9   Frank
           10  Gaussian
           11  Gumbel
           19  t


## Usage
               tau = CopulaKendallTau(family,theta)


## Inputs
       family      = The copula family.
       theta       = The parameter vector for the pair-copula.


## Outputs
      tau          = The value of Kendall's tau.
