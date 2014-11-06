# PairCopulaVfun

Computing the v-function of a pair-copula

---
## Purpose
The function computes the h-function, i.e., it evaluates the
        conditional distribution function of U2|U1 at u2, given U1=u1.
        Possible copula families:

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
               U = PairCopulaVfun(family,u1,u2,theta)
           Rotated pair-copulas
               U = PairCopulaVfun(family,u1,u2,theta,rotation)


## Inputs
       family    = The copula family.
       u1        = A (n x 1) dimensional vector of values lying in [0,1].
       u2        = A (n x 1) dimensional vector of values lying in [0,1].
       theta     = The parameter vector for the pair-copula.
       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
                   rotation is achieved by letting the rotation argument
                   empty or by choosing 0 rotation.


## Outputs
      U          = n dimensional vector of probabilities (sometimes 
                   called pseudo-observations), which are lying strictly
                   in [0,1].
