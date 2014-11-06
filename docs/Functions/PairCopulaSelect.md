# PairCopulaSelect

Selecting pair-copula families

---
## Purpose
The function selects pair-copula families by applying the test for
        independence of copula data of Genest and Favre (2007) using a
        significance level of 5  and afterwards, in the case of an
        rejection of the H0, it selects the "best" fitting pair-copula by
        optimizing with respect to Akaike's information criterion (AIC).
        The set of pair-copula families that are tested can be specified
        in several ways. By chosing the option 'all', the possible
        pair-copulas are:
        
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
           
Alternatively, one can specify a own subset of copula families or
        one can choose the option 'R-package', which coincides with the
        pair-copula families set of the R-package VineCopula of 
        Schepsmeier, Stöber, and Brechmann (2013). The possible rotation
        degrees for the rotatable pair-copula families are always choosen
        according to the empirical Kendall's tau values.


## Usage
           Pair copula selection from the whole set of possible families
           (default).
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2)
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2,'all')
           Pair copula selection from a user specified set of possible
           pair-copula families.
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2,familyset)
           Pair copula selection from the set of pair-copula families,
           which coincides with the one of the  R-package VineCopula of
           (cf. Schepsmeier, Stöber, and Brechmann (2013)), where the
           following commands are all equivalent.
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2,'R')
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2,'R-package')
       [family,ParamHat,rotation] = PairCopulaSelect(u1,u2,'VineCopulaPackage')


## Inputs
       u1        = A (n x 1) dimensional vector of values lying in [0,1].
       u2        = A (n x 1) dimensional vector of values lying in [0,1].
       familyset = The set of possible pair-copula families. By setting it
                   to the strings 'all', 'R', 'R-package' or
                   'VineCopulaPackage' one can choose one of the pre-
                   defined sets. Alternatively, one can choose an vector
                   containing a subset of the possible families.


## Outputs
       family    = The selected pair-copula family.
       ParamHat  = The maximum-likelihood estimates for the parameters of
                   the selected pair-copula family.
       rotation  = The selected rotation level (can of course be empty).


## References
Schepsmeier, U., J. Stöber and E. C. Brechmann (2013), VineCopula:
      Statistical inference of vine copulas, R package version 1.2, url:
      http://CRAN.R-project.org/package=VineCopula.
