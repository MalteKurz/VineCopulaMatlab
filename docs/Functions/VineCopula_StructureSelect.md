# StructureSelect (VineCopulaObject)

Selecting the structure and pair-copula families for a vine copula

---
## Purpose
The function can be used to find an adequate structure
        and to select pair-copula families for a given data
        set. By default (StructuringRule = 0) the nodes of
        the C-vine are chosen in a way that in each tree,
        the root (i.e. the node, which is connected by a
        copula to all other nodes) is the variable, which has
        maximal dependence with all other variables. The
        maximal dependence is found by choosing the variable,
        which has the maximal column sum in the matrix of
        absolute empirical Kendall’s τ (cf. Schepsmeier,
        Stöber, and Brechmann (2013) for an R-function
        (RVineStructureSelect) of exactly the same procedure
        and Czado, Schepsmeier, and Min (2012, p. 240) for the
        theoretical background of the approach).
        
Alternative methods for choosing a structure for
        C-Vine copulas (cf. Nikoloulopoulos, Joe and Li
        (2012)):
        The root node of the first tree is chosen in same way
        as in the default method, i.e., the variable with the
        strongest dependence with all other variables. Then
        one can choose between three rules suggested in
        Nikoloulopoulos, Joe and Li (2012, p.3665):
        StructuringRule = 1: List the other variables by their
                             dependence to the root node in
                             decreasing order.
        StructuringRule = 2: List the other variables by their
                             dependence to the root node in
                             increasing order.
        StructuringRule = 3: List the other variables
                             sequentially by choosing the
                             variable which is least dependent
                             with the previously selected one.

Furthermore, the copula families are chosen
        according to the AIC criterion and for each pair-
        copula an independence test is performed (cf.
        Schepsmeier, Stöber, and Brechmann (2013) and
        Brechmann and Schepsmeier (2013) for R-functions
        (RVineStructureSelect / RVineCopSelect /
        CDVineCopSelect) and Genest and Favre (2007, p. 351)
        for the independence test).

For D-Vine copulas there is no structuring rule
        implemented yet. Therefore, the function uses the
        specified structure and only selects the pair-copula
        families using the AIC criterion.


## Usage
        Select a simplified vine copula model
          VineCopulaHat = StructureSelect(VineCopulaObject,u)
        Select a simplified vine copula model, where the
        pair-copulas are from a specified set of families
        (possible choices 'all' (default), 'R', 'R-package',
        'VineCopulaPackage' (they all are equivalent and
        correspond to the set of pair-copulas of the R-package
        VineCopula (Schepsmeier, Stöber, and Brechmann (2013),
        Version 1.2)) or a cell-array consisting of possible
        pair-copula families (i.e., a user selected list of
        possible pair-copula families).
          VineCopulaHat = StructureSelect(VineCopulaObject,u,'all')
          VineCopulaHat = StructureSelect(VineCopulaObject,u,'R')
          VineCopulaHat = StructureSelect(VineCopulaObject,u,'R-package')
          VineCopulaHat = StructureSelect(VineCopulaObject,u,'VineCopulaPackage')
          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset)
        Select a simplified vine copula model using an
        alternative structuring rule:
          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,1)
          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,2)
          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,3)



## Inputs
        VineCopulaObject= An object from the class VineCopula.
        u               = A (n x d) dimensional vector of
                          values lying in [0,1] (the
                          observations).
        familyset       = The set of possible pair-copula
                          families. By setting it to the
                          strings 'all', 'R', 'R-package' or
                          'VineCopulaPackage' one can choose
                          one of the pre-defined sets.
                          Alternatively one can choose an
                          array containing a subset of the
                          possible families:
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


## Outputs
        VineCopulaHat   = An object from the class VineCopula.
                          The select vine copula structure can
                          be found in VineCopulaHat.structure
                          and the selected pair-copulas in
                          VineCopulaHat.families. Furthermore,
                          the sequential estimates, which are
                          obtained during the selection
                          procedure of the structure and
                          copula families are stored in
                          VineCopulaHat.parameters.


## References
 [1]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling
      Dependence with C- and D-Vine Copulas: The R-Package
      CDVine", Journal of Statistical Software 52(3), R
      package version 1.1-13, pp. 1-27, url:
      http://CRAN.R-project.org/package=CDVine.
      
 [2]  Czado, C., U. Schepsmeier, and A. Min (2012), "Maximum
      likelihood estimation of mixed C-vines with application
      to exchange rates", Statistical Modelling 12(3), pp.
      229-255.
      
 [3]  Genest, C. and A. Favre (2007), "Everything You Always
      Wanted to Know about Copula Modeling but Were Afraid to
      Ask", Journal of Hydrologic Engineering 12(4), pp. 347-
      368.
      
 [4]  Nikoloulopoulos, A. K., H. Joe, H. Li (2012), "Vine
      copulas with asymmetric tail dependence and applications
      to financial return data", Computational Statistics &
      Data Analysis 56(11), pp. 3659-3673.
      
 [5]  Schepsmeier, U., J. Stöber, and E. C. Brechmann (2013),
      VineCopula: Statistical inference of vine copulas, R
      package version 1.2, url:
      http://CRAN.R-project.org/package=VineCopula.
