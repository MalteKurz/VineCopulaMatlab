# VineCopulaRand

Generating (pseudo-)random variables from a copula or vine copula

---
## Purpose
The function draws N pseudo-random d-dimensional tuples from a
        C-Vine or D-Vine copula. The vine copulas can be of the simplified
        and non-simplified form. In the case of a non-simplified Vine
        copula, the parameters of the conditional copulas have to be given
        as functionals of the variables on which the conditioning is done.


## Usage
           Standard C-Vine copula (simplified & non-simplified, where the
           parameters of the conditional copulas have to be given as
           functionals of the variables on which the conditioning is done)
             simplified
               U = CopulaRand('C-Vine',N,d,families,thetas,ones(1,(d-1)*(d-2)/2))
             non-simplified
               U = CopulaRand('C-Vine',N,d,families,thetas,simplified,condparameterfunctionals)
           Standard D-Vine copula (simplified & non-simplified, where the
           parameters of the conditional copulas have to be given as
           functionals of the variables on which the conditioning is done)
             simplified
               U = CopulaRand('D-Vine',N,d,families,thetas,ones(1,(d-1)*(d-2)/2))
             non-simplified
               U = CopulaRand('D-Vine',N,d,families,thetas,simplified,condparameterfunctionals)


## Inputs
       type      = The vine copula type.
       N         = The number of pseudo-random tuples to be drawn.
       d         = The dimension of the C- or D-Vine.
       rotation    = A vector of the same dimension as families in
                     which one can specify rotation levels.
       families  = A cell array of the pair-copula families, which are
                   part of the PCC. The cell array has to have the length
                   (d-1)*d/2. The first d-1 entries are the copula
                   families in the first tree and the next d-2 entries are
                   the copula families in the second tree and so on. That
                   means, for d=4 the array should look similar to this
                   {'Frank', 'Frank', 'Frank', 'AMH', 'AMH', 'Clayton'},
                   which is the special case where all copulas in the
                   second tree are AMH copulas and all copulas in the
                   third tree are Clayton copulas. The order of the
                   families is (exemplarily for the case d=4): C12, C13,
                   C14, C23|1, C24|1, C34|12.
                   Note: If families is a simple string/character, e.g.,
                   'Clayton', then all pair-copulas are specified to be
                   from this copula family.
       thetas    = The values of the parameters for the (d-1)*d/2 pair-
                   copulas. These parameter values have to be given in the
                   same order as the families vector, but in a row-
                   vector. If a pair-copula is a independence copula, then
                   there is no parameter needed. Furthermore, if a pair-
                   copula has to or more parameters, the parameters have
                   to be given in same order as they have to be provided
                   if the pair-copula is considered only. For example for
                   a t-copula, the first parameter is rho and the second
                   parameter is the degrees of freedom parameter nu.
                   Note: If thetas is a skalar then all parameters (i.e.,
                   the parameters of all pair-copulas) are set to the same
                   value (i.e. to the one given).
       simplified
                 = A vector consisting of (d-1)*(d-2)/2 zeros or ones
                   (either 1 (for an unconditional bivariate pair-copula)
                   or 0 for a conditional copula) specifying, whether
                   the copulas in the second and higher tree are pair-copulas
                   (unconditional bivariate copulas) or conditional copulas.
       condparameterfunctionals
                 = A cell array consisting of (# non-partial (i.e.
                   conditional) pair-copulas) cell array, which consist of
                   as many parameter functionals as the copula families
                   has parameters. Therefore, if one of the parameters is
                   not varying in the conditioning set, then the function
                   handle has to give back constant (e.g., @(w) c, where c
                   is a scalar constant. Functionals for conditional
                   copulas beeing part of the second tree of the C-Vine
                   have to be functions of one variable (i.e., they should
                   be evaluatable for column vectors). Functionals of the
                   third tree of the C-Vine have to be functionals of two
                   variables (i.e. they should be evaluatable for matrices
                   of the size N x 2). And so on, for conditional copulas,
                   where the conditioning is done on three or more
                   variables. The ordering of the variables in the
                   conditioning set corresponds to the ordering of the
                   root nodes of the different trees (e.g., in the fourth
                   tree, the first conditional copula is the C_45|123
                   copula and the conditioning set is three-dimensional,
                   where the first column of observations corresponds to
                   observations of the root node in the first tree (i.e.
                   variable X_1).
                   Note: Of course it is possible, that the parameter of a
                   conditional copula only depends on some of the
                   variables, but then the functions have to be specified
                   in an appropriate way (e.g., if a conditional copula in
                   the third tree, where the conditioning is done on two
                   variables, only depends on the second variable, the
                   function could look similar to this @(u) 3.*u(:,2),
                   when the parameter of the conditional copula is exactly
                   three times the second variable on which the
                   conditioning is done.


## Outputs
      U          = (N x d) Matrix of simulated tuples from the specified
                   copula, where every row is a simulated d-dimensional
                   tuple from the d-dimensional C-Vine or D-Vine copula.
