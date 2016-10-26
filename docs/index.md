# VineCopulaMatlab toolbox
___________________________________________________________________________
A **MATLAB toolbox** for vine copulas based on the C++ shared library VineCopulaCPP
___________________________________________________________________________

## Description of the Vine Copulas with C++ toolbox
The toolbox can be used for high-dimensional dependence modeling with vine copula models. A key feature of the toolbox is a framework, which allows to test whether the simplifying assumption is a reasonable assumption for approximating high-dimensional distributions using simplified vine copula models.

## Highlights
  * Modeling of (high-)dimensional [0,1]-data by C-vine and D-vine copulas.
  * 20 different pair-copula families (62 families with rotated pair-copulas).
  * The most important object class VineCopula is implemented in MATLAB.
  * Functions for simulating from simplified and non-simplified C- and D-vine copulas. In the case of the non-simplified C- or D-vines, the parameters of all conditional bivariate copulas can be specified as functions of the conditioning variables.
  * Functions for selecting and (jointly-)estimating C- and D-vine copula models.
  * Functions for (jointly-)estimating C- and D-vine copula models.
  * A vectorial independence test, which can be used for testing sequentially the simplified assumption for C- and D-vine copulas.
  * Most computations are implemented in C++. Parts can be optionally performed parallel.

## Hosting and bug reporting
  * The VineCopulaMatlab toolbox is hosted at GitHub and can be found under [https://github.com/MalteKurz/VineCopulaMatlab](https://github.com/MalteKurz/VineCopulaMatlab).
  * The C++ shared library VineCopulaCPP is also hosted at GitHub and can be found under [https://github.com/MalteKurz/VineCopulaCPP](https://github.com/MalteKurz/VineCopulaCPP).
  * Please use GitHub issues ([VineCopulaMatlab issues](https://github.com/MalteKurz/VineCopulaMatlab/issues) and [VineCopulaCPP issues](https://github.com/MalteKurz/VineCopulaCPP/issues) for reporting any issues or bugs.

## Demo
Please see the [demo](DemoVineCPP.md) for further details about the functionality of the VineCPP toolbox.

## Remarks
  * The central class of the toolbox is the VineCopula class. Working with this class, for most input variables consistency checks are performed.
  * In all the other functions (primarily the functions for two-dimensional
    (pair-)copulas there are not that many consistency checks performed.
  * At the moment within this toolbox you can only choose between C-vines and D-vines. The superclass of regular vines (R-vines) is not implemented yet. If you also want to use R-Vines then the R-package VineCopula (cf. Schepsmeier, Stöber, and Brechmann (2013)) is an excellent alternative. 
 
## Dependencies
  * A **C++-compiler** being compatible with the used MATLAB release (cf.
    http://www.mathworks.de/support/compilers for informations about
    compatible compilers for different releases and opperating systems).
  * The C++ libraries **boost** (cf. http://www.boost.org/). Used are some functions for statistical distributions from the Boost Math Toolkit. Furthermore, the boost libraries are used for random number generation.
  * The nonlinear optimization library **NLopt** (http://ab-initio.mit.edu/wiki/index.php/NLopt).
  * **OpenMP** for parallel computing (http://openmp.org/wp/).
  * The Fortran 77 routine **MVTDST** (file mvtdstpack.f) from (http://www.math.wsu.edu/faculty/genz/software/software.html; Alan Genz). (It is only needed for computing the CDF of the bivariate normal and t copula.)
  * The C++ library **VineCopulaCPP** (https://github.com/MalteKurz/VineCopulaCPP).
 
## References
 [1]  Aas, K., C. Czado, A. Frigessi and H. Bakken (2009), "Pair-copula constructions of multiple dependence", Insurance: Mathematics and Economics 44(2), pp. 182-198.

 [2]  Acar, E. F., C. Genest and J. Neslehová (2012), "Beyond simplified pair-copula constructions", Journal of Multivariate Analysis 110, pp. 74-90.

 [3]  Balakrishnan, N. and Lai, C.-D. (2009), "Continuous Bivariate Distributions", 2. ed., New York, NY: Springer.

 [4]  Bedford, T. and R. M. Cooke (2001), "Probability density decomposition for conditionally dependent random variables modeled by vines", Annals of Mathematics and Artificial Intelligence 32 (1), pp. 245-268.

 [5]  Bedford, T. and R. M. Cooke (2002), "Vines -- A new graphical model for dependent random variables", The Annals of Statistics 30(4), pp. 1031-1068.

 [6]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling Dependence with C- and D-Vine Copulas: The R-Package CDVine", Journal of Statistical Software 52(3), R package version 1.1-13, pp. 1-27, url: http://CRAN.R-project.org/package=CDVine.

 [7]  Czado, C., U. Schepsmeier, and A. Min (2012), "Maximum likelihood estimation of mixed C-vines with application to exchange rates", Statistical Modelling 12(3), pp. 229-255.

 [8]  Eschenburg, P. (2013), "Properties of extreme-value copulas", Diploma Thesis, Fakultät für Mathematik, Technische Universität München, url: http://mediatum.ub.tum.de/download/1145695/1145695.pdf.

 [9]  Genest, C. and A. Favre (2007), "Everything You Always Wanted to Know about Copula Modeling but Were Afraid to Ask", Journal of Hydrologic Engineering 12(4), pp. 347-368.

 [10] Hobaek-Haff, I., K. Aas and A. Frigessi (2010), "On the simplified pair-copula construction -- Simply useful or too simplistic?", Journal of Multivariate Analysis 101(5), pp. 1296-1310.

 [11] Joe, H. (1996), "Families of m-Variate Distributions With Given Margins and m(m-1)/2 Bivariate Dependence Parameters", Distributions with Fixed Marginals and Related Topics, ed. by L. Rüschendorf, B. Schweizer, and M. D. Taylor, Hayward, CA: Institute of Mathematical Statistics.

 [12] Joe, H. (1997), Multivariate models and dependence concepts, 1. ed., reprint., Monographs on statistics and applied probability; 73, Boca Raton, Fla. [u.a.]: Chapman & Hall/CRC.

 [13] Kojadinovic, I. and M. Holmes (2009), "Tests of independence among continuous random vectors based on Cramér-von Mises functionals of the empirical copula process", Journal of Multivariate Analysis 100(6), pp. 1137-1154.

 [14] Kosorok, M. R. (2008), Introduction to Empirical Processes and Semiparametric Inference, Springer Series in Statistics, New York, NY: Springer.

 [15] Kurowicka, D. and H. Joe (Eds.) (2011), "Dependece Modeling -- Vine Copula Handbook", Singapore: World Scientific Publishing Co. Pte. Ltd.

 [16] Nelsen, R. B. (2006), "An introduction to copulas", 2. ed., Springer series instatistics, New York, NY: Springer.

 [17] Omelka, M. and M. Pauly (2012), "Testing equality of correlation coefficients in two populations via permutation methods, Journal of Statistical Planning and Inference 142, pp. 1396-1406.

 [18] Patton, A. J. (2002), "Applications of Copula Theory in Financial Econometrics", Unpublished Ph.D. dissertation, University of Colifornia, San Diego, url: http://www.amstat.org/sections/bus_econ/papers/patton_dissertation.pdf.

 [19] Patton, A. J. (2006), "Modelling asymmetric exchange rate dependence", International Economic Review 47(2), pp. 527-556.

 [20] Quessy, J.-F. (2010), "Applications and asymptotic power of marginal-free tests of stochastic vectorial independence", Journal of Statistical Planning and Inference 140(11), pp. 3058-3075.

 [21] Rémillard, B. (2013), "Statistical Methods For Financial Engineering", Boca Raton, FL: Chapman & Hall.

 [22] Rémillard, B. and O. Scaillet (2009), "Testing for equality between two copulas", Journal of Multivariate Analysis 100(3), pp. 377-386.

 [23] Schepsmeier, U., J. Stöber and E. C. Brechmann (2013), VineCopula: Statistical inference of vine copulas, R package version 1.2, url: http://CRAN.R-project.org/package=VineCopula.

 [24] Segers, J. (2012), "Asymptotics of empirical copula processes under non-restrictive smoothness assumptions", Bernoulli 18(3), pp. 764-782.
 
 [25] Spanhel, F., Kurz, M.S., 2015. Simplified vine copula models: Approximations based on the simplifying assumption. ArXiv e-prints https://arxiv.org/abs/1510.06971.

 [26] Stöber, J., H. Joe and C. Czado (2013), "Simplified pair copula constructions -- Limitations and extensions", Journal of Multivariate Analysis 119, pp. 101-118.

 [27] van der Vaart, A. W. and J. A. Wellner (1996), Weak Convergence and Empirical Processes -- With Applications to Statistics, Springer Series in Statistics, New York [u.a.]: Springer.


___________________________________________________________________________
 Author: Malte Kurz
 
___________________________________________________________________________
