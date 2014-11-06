# PairCopulaRand

Generating (pseudo-)random variables from a pair-copula

---
## Purpose
The function draws N pseudo-random 2-dimensional tuples from a
        pair-copula. Possible copula families:
          
           0   Indep
           1   AMH
           2   AsymFGM
           3   BB1
           4   BB6
           5   BB7
           6   BB8
           7   Clayton
           8   FGM
           9   Frank
           10  Gaussian
           11  Gumbel
           12  IteratedFGM
           13  Joe
           14  PartialFrank
           15  Plackett
           16  Tawn1
           17  Tawn2
           18  Tawn
           19  t


## Usage
               U = CopulaRand(family,N,theta)
           Rotated pair-copulas
               U = CopulaRand(family,N,theta,rotation)


## Inputs
       family    = The copula family.
       N         = The number of pseudo-random tuples to be drawn.
       theta     = The parameter vector for the pair-copula.
       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
                   rotation is achieved by letting the rotation argument
                   empty or by choosing 0 rotation.


## Outputs
      U          = (N x 2) Matrix of simulated tuples from the specified
                   pair-copula.
