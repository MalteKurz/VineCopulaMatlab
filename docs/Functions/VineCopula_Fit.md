# Fit (VineCopulaObject)

Estimating objects of the VineCopula class

---
## Purpose
The function computes ML-estimates for the parameters
        of a simplified vine copula. Therefore, first starting
        values for the joint estimation are obtained by
        iteratively estimating the pair-copulas in the first
        trees and using those estimates to obtain the
        arguments for the copulas in the second tree. Then the
        pair-copulas in the second tree are estimated and so
        on. These estimated parameters from the sequential
        procedure are then used to obtained the ML-estimates,
        by minimizing the overall negative log-likelihood of
        the whole vine copula numerically.


## Usage
        Estimating a simplified vine copula (joint estimation;
        the default method)
          VineCopulaHat = Fit(VineCopulaObject,u)
          VineCopulaHat = Fit(VineCopulaObject,u,'joint')
        Estimating a simplified vine copula (sequential
        estimation)
          VineCopulaHat = Fit(VineCopulaObject,u,'sequential')
        Estimating a simplified vine copula (with a cut off
        tree / truncation level)
          VineCopulaHat = Fit(VineCopulaObject,u,EstMethod,CutOffTree)


## Inputs
        VineCopulaObject= An object from the class VineCopula.
        u               = A (n x d) dimensional vector of
                          values lying in [0,1] (the
                          observations).
        EstMethod       = The estimation method must be either
                          'joint' or 'sequential'. If it is
                          not explicitly given, a joint
                          estimation is performed (default).
        CutOffTree      = The CutOffTree (or also called
                          truncation level) can be used to set
                          all pair-copulas from the (CutOffTree
                          + 1)-th tree on to independence
                          copulas (i.e., ignore them in the
                          joint estimation). The CutOffTree
                          does only influence the joint
                          estimation.


## Outputs
       VineCopulaHat   = An object from the class VineCopula.
                         The sequential estimates are stored
                         in VineCopulaHat.SeqEstParameters,
                         the estimated parameters from the
                         joint estimation are stored in
                         VineCopulaHat.parameters and the two
                         maximized values of the vine copula
                         log-likelihood are stored in
                         VineCopulaHat.MaxLLs.
