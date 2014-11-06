
---

# Installation

For using the VineCPP toolbox you should first check whether all dependencies are installed.

## Dependencies
  * A **C++-compiler** being compatible with the used MATLAB release (cf.
    http://www.mathworks.de/support/compilers for informations about
    compatible compilers for different releases and opperating systems).
  * The C++ libraries **boost** (cf. http://www.boost.org/). Used are some functions for statistical distributions from the Boost Math Toolkit. Furthermore, the boost libraries are used for random number generation.
  * The nonlinear optimization library **NLopt** (http://ab-initio.mit.edu/wiki/index.php/NLopt).
  * **OpenMP** for parallel computing (http://openmp.org/wp/).
  * The Fortran 77 routine **MVTDST** (file mvtdstpack.f) from (http://www.math.wsu.edu/faculty/genz/software/software.html; Alan Genz). (It is only needed for computing the CDF of the bivariate normal and t copula.)
  
With the following code you can check whether you have installed the needed software (C++ compiler, boost libraries, NLopt and OpenMP) on your UNIX system.

    $ dpkg -s g++ | grep 'Version'
    4:4.8.2-1ubuntu6
    $ dpkg -s libboost-all-dev | grep 'Version'
    1.54.0.1ubuntu1
    $ dpkg -s libnlopt-dev | grep 'Version'
    2.4.1+dfsg-1ubuntu1
    dpkg -s libgomp1 | grep 'Version'
    4.8.2-19ubuntu1

## Getting started / Setting up the VineCPP toolbox
 1. Check whether you have a compatible **C++ compiler** (cf. http://www.mathworks.de/support/compilers for information about compatible compilers for different releases and operating systems). By the command mex -setup the compiler can be set.
 2. Install the C++ libraries **boost** (cf. http://www.boost.org/).
 3. Install the nonlinear optimization library **NLopt** (http://ab-initio.mit.edu/wiki/index.php/NLopt).
 4. Optionally parallel computing with **OpenMP** can be used. Note that your C++ compiler needs to support parallel computing with OpenMP.
 5. Download the Fortran 77 routine **MVTDST** (file mvtdstpack.f) from (http://www.math.wsu.edu/faculty/genz/software/software.html; Alan Genz). (It is only needed for computing the CDF of the bivariate normal and t copula.)
 6. Install the VineCPP toolbox by executing the installVineCPP function.
  
## Remarks  
  * Especially if you are not working with an UNIX operating system you may
    get errors when compiling the C++ files during the execution of the
    installVineCPP function. You may have to add the path of the installation
    of the boost libraries by changing line 63 of the installVineCPP function
    from eval(['mex ' CPP_files{i}]) to eval(['mex -I''pathname'' ' CPP_files{i}]) , where pathname is replaced by the path of the boost libraries.
  * For some C++ compilers the function expm1 is not available. A
    possible workaround is to use the expm1 from the boost libraries.
    Therefore, add the line #include <boost/math/special_functions/expm1.hpp> 
    to the header file VineCPP_header.hpp. Then substitute expm1 by boost::math::expm1.
  * The function lgamma is also not available in some C++ compilers. A
    possible workaround is to use the lgamma function from the boost libraries.
    Therefore, add the line #include <boost/math/special_functions/gamma.hpp> 
    to the header file VineCPP_header.hpp. Then substitute lgamma
    by boost::math::lgamma.
    
---
    
# Handling of seeds for random number generation
**Boost.Random** is used for random number generation in C++ (cf. http://www.boost.org/doc/libs/1_56_0/doc/html/boost_random.html). If you execute installVineCPP there will be generated some seed for random number generation. With the function SetSeedVineCPP you can always set a new one. Further note:

* The seed is saved in the file "/private/Seed.dat".
* The seed is generate using the current system time.
* If you want to save the seed for future usage, please make a copy
          of the file "Seed.dat". Every time random numbers are generated in
          C++ the file will be overwritten.
* If the Seed should be set back to some former state, you just have
          to replace the "Seed.dat" file in "/private/" by your copy of
          Seed.dat you have previously generated.
* For further information please visit:
           http://www.boost.org/doc/libs/1_56_0/doc/html/boost_random/reference.html

---

# Parallel computing within the VineCPP toolbox using OpenMP
The following functions can be (partly) computed in parallel:

* PairCopulaSelect
* StructureSelect(VineCopula)
* Fit(VineCopula)
* VineCopulaFit
* VineCopulaNegLL

For effective usage of the parallelization you should enable nested parallel regions by:
    
    export OMP_NESTED=TRUE

The maximum number (N) of threads generated can be set via:

    export OMP_NESTED=TRUE
  
---