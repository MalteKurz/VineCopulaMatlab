# PairCopulaAIC

Computing the AIC of a pair-copula

---

##Purpose
The function computes the value of the AIC for a copula for a
        given matrix of observations u, which have to lie in the
        2-dimensional unit cube, evaluated at the ML estimates. Possible
        pair-copula families:
           
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
        [AIC,ParamHat] = PairCopulaAIC(family,u1,u2)
      Rotated pair-copulas
        [AIC,ParamHat] = PairCopulaAIC(family,u1,u2,roatation)


## Inputs
       family          = The copula family.
       u1              = A (n x 1) dimensional vector of values lying in [0,1].
       u2              = A (n x 1) dimensional vector of values lying in [0,1].
       rotation        = The degree of rotation, i.e., either 90, 180 or
                         270. No rotation is achieved by letting the
                         rotation argument empty or by choosing 0
                         rotation.


## Outputs
      AIC              = The value of the AIC evaluated at the ML-estimator.
      ParamHat         = The ML-estimate for the parameter vector.
