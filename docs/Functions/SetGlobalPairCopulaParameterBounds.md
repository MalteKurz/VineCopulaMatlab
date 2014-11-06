# SetGlobalPairCopulaParameterBounds

Setting the global (i.e., for the whole VineCPP toolbox) pair-copula parameter bounds, which are used for estimation and consistency checks

---
# Purpose
The GUI can be used to set the global (i.e., for the whole VineCPP
        toolbox) pair-copula parameter bounds, which are used for
        estimation and consistency checks. Possible copula families (the
        default values for the parameter bounds are given in parenthesis,
        where the ordering is as follows (lb1, ub1; lb2, ub2; lb3, ub3)):

           0   Indep
           1   AMH (-1, 1)
           2   AsymFGM (0, 1)
           3   BB1 (0, 6; 1, 6)
           4   BB6 (1, 6; 1, 6)
           5   BB7 (1, 6; 0.001, 6)
           6   BB8 (1, 6; 0.001, 1)
           7   Clayton (0, Inf)
           8   FGM (-1, 1)
           9   Frank (-30, 30)
           10  Gaussian (-0.999, 0.999)
           11  Gumbel (1, Inf)
           12  IteratedFGM (-1, 1; -1, 1)
           13  Joe (1, Inf)
           14  PartialFrank (0, 30)
           15  Plackett (0.001, Inf)
           16  Tawn1 (1.001, 20; 0.001, 0.999)
           17  Tawn2 (1.001, 20; 0.001, 0.999)
           18  Tawn (1.0001, 20; 0.001, 0.999; 0.001, 0.999)
           19  t (-0.999, 0.999; 1, 30)


# Usage
               SetGlobalPairCopulaParameterBounds
