# PairCopulaHVolume

Computing the h-volume of a copula

---
## Purpose
The function computes the probability for a random vector, being
        distributed according to a specific copula, to lie in a
        hyperrectangle. The hyperrectangle is defined by the cartesean
        product of the intervals specified by the lower bounds a and upper
        bounds b. Possible pair-copula families:
        
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
               P = CopulaHVolume(family,a,b,theta)
           Rotated pair-copulas
               P = CopulaHVolume(family,a,b,theta,rotation)


## Inputs
       family    = The copula family.
       a         = A 2-dimensional vector of lower bounds for 2 intervals,
                   defining a 2-dimensional hyperrectangle.
       b         = A 2-dimensional vector of upper bounds for 2 intervals,
                   defining a 2-dimensional hyperrectangle.
       theta     = The parameter vector for the pair-copula.
       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
                   rotation is achieved by letting the rotation argument
                   empty or by choosing 0 rotation.


## Outputs
      P          = The probability for a random variable, which is
                   distributed according to the specified pair-copula, to
                   lie in the hyperrectangle defined by the lower and
                   upper bounds vectors a and b.
