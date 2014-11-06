# VineCopula class

The class of vine copula models (yet C-Vines and D-Vines only)

---
## Properties
### dimension
The dimension (=d) of the vine copula.
###   type
The type of the vine copula, i.e.,
                              'C-Vine', 'CVine' or 0 and 'D-Vine'
                              ,'DVine' or 1.
###   simplified
A logical or vector of logicals, which
                              provides the information, whether
                              the simplifying assumption is
                              fulfilled or not. Either the whole
                              vine is set to simplified, i.e., all
                              conditional copulas are chosen to be
                              unconditional bivariate copulas (pair-
                              copulas) or each conditional copula
                              being part of the vine copula
                              model is set to 1 for a pair-copula,
                              i.e., a unconditional bivariate
                              copula or 0 for a conditional
                              copula.
                              
Possible values for simplified are 1
                              (for unconditional bivariate
                              copulas; pair-copulas) and 0 for
                              conditional copulas. If simplified
                              is set to 1 then every pair-copula
                              being part of the PCC is set to be a
                              partial copula. In contrast, if the
                              simplified property is chosen to be
                              0 then every conditional copula is
                              set to be a conditional copula with
                              functional parameter (note that this
                              can also be a unconditional copulas
                              whenever the functional parameter is
                              a constant). Alternatively, the
                              simplified property can also be
                              given as a 1 x ((d-1)*(d-2)/2)
                              matrix, where each entry is either
                              0 or 1 and therefore every copula is
                              set to be a pair-copula or a
                              conditional copula with functional
                              parameter, separately.
###   structure
A vector containing the information
                              about the structure of the vine
                              copula model.
                              
Possible values for the structure in
                              the case of a C-Vine copula are all
                              permutations of the (1,...,d)
                              sequence. The first entry of the
                              structure is the root node of tree 1
                              (i.e., the unique node of degree d-1),
                              the second entry is the root node of
                              degree d-2 in tree 2, etc.
                              For a D-Vine possible values for the
                              structure vector are also all
                              permutations of the (1,...,d)
                              sequence. The choice [1 2 3 4 5]
                              corresponds to the D-Vine consisting
                              of the pair-copulas C_12, C_23,
                              C_34, C_45, C_13|2, C_24|3, C_35|4,
                              C_14|23, C_25|34 and C_15|234.
### families
A 2 x (d*(d-1)/2) cell array or
                              matrix, where in the first row the
                              copula family for each pair-copula
                              is stored and in the second row, one
                              can optionally store the degree of
                              rotation for the pair-copula.
                              
Possible values for families is either
                              a 1 x (d*(d-1)/2) cell array (if there
                              are no rotated pair-copulas within the
                              PCC) or a 2 x (d*(d-1)/2) cell array,
                              where the first row consists of the
                              pair-copula families and the second
                              one of the degrees of rotation of the
                              different pair-copulas. If rotation is
                              empty for a pair-copula, then it is
                              the standard (unrotated) copula of
                              the specified family.
                              The copula families and rotation
                              degrees can also be given in the
                              numerical coding (e.g. 7 for the
                              Clayton copula and 90 for 'r90').
                              The order of the entries into the
                              columns is tree by tree and within
                              the tree it is ordered according to
                              the structure. For example if the
                              structure is [2,1,3,4] then the
                              entries in the family cell array
                              correspond to the pair-copulas in the
                              following ordering (for the C-Vine):
                              C_21, C_23, C_24, C_13|2, C_14|2 and
                              C_34_21. Here it is important to note
                              that the root in each tree is always
                              the first variable of each pair-
                              copula, which is important for the
                              specification of pair-copulas with
                              asymetric dependence. For a
                              five-dimensional D-Vine copula with
                              structure [1 2 3 4 5], the ordering
                              of the families has to correspond to
                              the pair-copulas C_12, C_23, C_34,
                              C_45, C_13|2, C_24|3, C_35|4,
                              C_14|23, C_25|34 and C_15|234.
###   parameters
A 1 x (d*(d-1)/2) cell array
                              containing all the parameters of the
                              unconditional pair-copulas.
                              Alternatively,  the same information
                              (the parameters) can also be given
                              simply in vectorized format, where
                              the parameters are contained in the
                              vector in the same ordering as for
                              example the pair-copula families.
                              
Possible values for parameters are
                              either a 1 x (d*(d-1)/2) cell array,
                              where each entry is a vector of the
                              corresponding pair-copula
                              parameter(s). Alternatively, one can
                              also specify the parameters by
                              vectorizing the information contained
                              in the cell array of parameters.
###   condparameterfunctionals
A 1 x ((d-1)*(d-2)/2) cell array
                              containing parameter functionals for
                              the conditional copulas being part
                              of the vine. The parameter
                              functionals can also be provided in
                              vectorized form.
                              
Possible values for the
                              condparameterfunctionals are either a
                              1 x ((d-1)*(d-2)/2) cell array, where
                              each entry is a vector of the
                              corresponding pair-copula parameter
                              functional(s). Alternatively one can
                              also specify the parameter
                              functional(s) by vectorizing the
                              information contained in the cell
                              array of parameter functional(s).
### MaxLLs
The values of the maximized log-
                              likelihoods of the PCC evaluated at
                              the sequentially and jointly
                              estimated parameter vectors.
                              
The first entry always corresponds to
                              the maximized value of the vine copula
                              log-likelihood evaluated at the
                              sequentially estimated parameters,
                              while the second entry corresponds to
                              the jointly estimated parameters.
### SeqEstParameters
The sequentially estimated
                              parameters for the vine copulas.
                              
The estimates from estimating the
                              vine copula by applying the
                              sequential estimation approach.

## Methods
###  [VineCopula (VineCopulaObject)](VineCopula_Constructor.md)
The "function" VineCopula can be
                              used to construct members of the
                              VineCopula class. At least the
                              properties dimension, type and
                              simplified have to be specified.
###   [Fit (VineCopulaObject)](VineCopula_Fit.md)
The method Fit can be used to
                              estimate a vine copula model,
                              specified as a object of the
                              VineCopula class, using a dataset by
                              maximizing the log-likelihood of the
                              vine copula.
###  [Sim (VineCopulaObject)](VineCopula_Sim.md)
The method Sim can be used to
                              simulate from a vine copula,
                              specified as a object of the
                              VineCopula class.
###  [StructureSelect (VineCopulaObject)](VineCopula_StructureSelect.md)
The method StructureSelect can be
                              used to find a "adequat" structure
                              and pair-copula families (and
                              parameters) for a given data set.
###  GetPseudoObsFromVine
The method GetPseudoObsFromVine can
                              be used to obtain pseudo-
                              observations from the conditional
                              copulas being part of the specified
                              PCC.
###  [SeqTestOnSimplified (VineCopulaObject)](VineCopula_SeqTestOnSimplified.md)
The method SeqTestOnSimplified can
                              be used to sequentially test on the
                              simplifying assumption by applying
                              a vectorial independence test.
    
## References
[1]  Aas, K., C. Czado, A. Frigessi and H. Bakken (2009), "Pair-copula constructions of multiple dependence", Insurance: Mathematics and Economics 44(2), pp. 182-198.

 [2]  Acar, E. F., C. Genest and J. Neslehová (2012), "Beyond simplified pair-copula constructions", Journal of Multivariate Analysis 110, pp. 74-90.
 
 [3]  Bedford, T. and R. M. Cooke (2001), "Probability density decomposition for conditionally dependent random variables modeled by vines", Annals of Mathematics and Artificial Intelligence 32 (1), pp. 245-268.
 
 [4]  Bedford, T. and R. M. Cooke (2002), "Vines -- A new graphical model for dependent random variables", The Annals of Statistics 30(4), pp. 1031-1068.
 
 [5]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling Dependence with C- and D-Vine Copulas: The R-Package CDVine", Journal of Statistical Software 52(3), R package version 1.1-13, pp. 1-27, url: http://CRAN.R-project.org/package=CDVine.
 
 [6]  Hobaek-Haff, I., K. Aas and A. Frigessi (2010), "On the simplified pair-copula construction -- Simply useful or too simplistic?", Journal of Multivariate Analysis 101(5), pp. 1296- 1310.
 
 [7]  Joe, H. (1996), "Families of m-Variate Distributions With Given Margins and m(m-1)/2 Bivariate Dependence Parameters", Distributions with Fixed Marginals and Related Topics, ed. by L. Rüschendorf, B. Schweizer, and M. D. Taylor, Hayward, CA: Institute of Mathematical Statistics.
 
 [8]  Joe, H. (1997), Multivariate models and dependence concepts, 1. ed., reprint., Monographs on statistics and applied probability; 73, Boca Raton, Fla. [u.a.]: Chapman & Hall/CRC.
 
 [9]  Kurowicka, D. and H. Joe (Eds.) (2011), "Dependece Modeling -- Vine Copula Handbook", Singapore: World Scientific Publishing Co. Pte. Ltd.
 
 [10] Kurz, M. (2012), "On the simplified pair-copula construction -- Graphical and computational illustrations", Unpublished Term Paper, Master Seminar on Financial Econometrics, Department of Statistics, Ludwig-Maximilians-University Munich.
 
 [11] Kurz, M. (2013), "Tests on the partial copula", Unpublished Master's Thesis, Department of Statistics, Ludwig-Maximilians-University Munich.
 
 [12] Nikoloulopoulos, A. K., H. Joe, H. Li (2012), "Vine copulas with asymmetric tail dependence and applications to financial return data", Computational Statistics & Data Analysis 56(11), pp. 3659-3673.
 
 [13] Schepsmeier, U., J. Stöber and E. C. Brechmann (2013), VineCopula: Statistical inference of vine copulas, R package version 1.2, url: http://CRAN.R-project.org/package=VineCopula.
 
 [14] Spanhel, F. and M. Kurz (2014), "Simplified vine copula approximations -- Properties and consequences", submitted for publication. 
 
 [15] Stöber, J., H. Joe and C. Czado (2013), "Simplified pair copula constructions -- Limitations and extensions", Journal of Multivariate Analysis 119, pp. 101-118.
