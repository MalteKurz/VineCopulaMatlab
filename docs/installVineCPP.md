
---

# Installation

For using the VineCopulaMatlab toolbox you should first check whether all dependencies are installed.

## Dependencies
  * A **C++-compiler** being compatible with the used MATLAB release (cf.
    http://www.mathworks.de/support/compilers for informations about
    compatible compilers for different releases and opperating systems).
  * The C++ libraries **boost** (cf. http://www.boost.org/). Used are some functions for statistical distributions from the Boost Math Toolkit. Furthermore, the boost libraries are used for random number generation.
  * The nonlinear optimization library **NLopt** (http://ab-initio.mit.edu/wiki/index.php/NLopt).
  * **OpenMP** for parallel computing (http://openmp.org/wp/).
  * The Fortran 77 routine **MVTDST** (file mvtdstpack.f) from (http://www.math.wsu.edu/faculty/genz/software/software.html; Alan Genz). (It is only needed for computing the CDF of the bivariate normal and t copula.)

---

# Installation under linux
With the following code you can check whether you have installed the needed software (C++ compiler, boost libraries, NLopt and OpenMP) on your UNIX system.

    $ dpkg -s g++ | grep 'Version'
    4:4.8.2-1ubuntu6
    $ dpkg -s libboost-all-dev | grep 'Version'
    1.54.0.1ubuntu1
    $ dpkg -s libnlopt-dev | grep 'Version'
    2.4.1+dfsg-1ubuntu1
    dpkg -s libgomp1 | grep 'Version'
    4.8.2-19ubuntu1

## Getting started / Setting up the VineCopulaMatlab toolbox
 1. Check whether you have a compatible **C++ compiler** (cf. http://www.mathworks.de/support/compilers for information about compatible compilers for different releases and operating systems). By the command mex -setup the compiler can be set.
 2. Install the C++ libraries **boost** (cf. http://www.boost.org/).
 
    ```$ apt-get install libboost-all-dev```
 
 3. Install the nonlinear optimization library **NLopt** (http://ab-initio.mit.edu/wiki/index.php/NLopt).
 
    ```$ apt-get install libnlopt-dev```
 
 4. Optionally parallel computing with **OpenMP** can be used. Note that your C++ compiler needs to support parallel computing with OpenMP.
 
    ```$ apt-get install libgomp1```
 
 5. Change the path to the directory you have saved the VineCopulaMatlab folder, i.e.,
 
    ```$ cd [path]/VineCopulaMatlab```
 
 6. Build the VineCopulaMatlab toolbox using the make command.
 
    ```$ make```
    
 7. Install the VineCopulaMatlab toolbox using the make command.
 
    ```$ sudo make install```
 
---
   
# Installation under Windows

## Setting up a C++ build environment for Windows
A nice step-by-step guidance can be found at [ascend4.org](http://ascend4.org/Setting_up_a_MinGW-w64_build_environment). Following the steps for setting up *MSYS*, *Switchable 32- and 64-bit modes* and *MinGW-w64* you will get a GCC (*TDM-GCC MinGW-w64*) for Windows including a command window. It also includes [OpenMP](http://openmp.org/wp/) for parallel computing.
## Getting started / Setting up the VineCopulaMatlab toolbox
 1. Download the `mex_C++_mingw-w64.xml` file from [Chappjc's GitHub folder MATLAB](https://github.com/chappjc/MATLAB/tree/master/MinGW). Open the file and add `-fopenmp -lgfortran -lgomp -lnlopt -lVineCopulaCPP ` at the end of line 67 (definition of the variable `LINKLIBS`) directly after `-lmex -lmx -leng -lmat -lmwlapack -lmwblas `. Save the file `mex_C++_mingw-w64.xml` in the `private` folder of the VineCopulaMatlab toolbox.
 2. The C++ libraries boost (cf. http://www.boost.org/) don't need to be build. Just add the include directory when building the toolbox.
 3. The [nonlinear optimization library NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt). If you have set up MingGW-w64 you can just follow the installation guidance, i.e.,
 
        ./configure --prefix=Installation_Dir
        make
        make install
        
 4. Build the VineCopulaMatlab toolbox using the make command. You may have to set the boostincludedir option. If the boost libraries files are in the folder `C:\BoostDir`, then you may set boostincludedir to /c/BoostDir.
 
    ```$ make boostincludedir=/c/BoostDir```
    
 5. Install the VineCopulaMatlab toolbox using the make command. Depending on the installation directory (set as prefix) you may need to run msys as admin. You may have to set the prefix and boostincludedir options. If you set prefix to `/c/Dir` the library files will be placed in `C:\Dir\lib` and the include files will be placed in `C:\Dir\include`. Additionally you may want to specify the variables `mingwroot` (default value = "C:\\MinGW\\64") and `mingw_mex_xml` (default value = "mex_C++_mingw-w64.xml"), which specify where you have installed MinGW-w64 and how your mex setup file is called.
 
    ```$ make install prefix=/c/Dir boostincludedir=/c/BoostDir```
         
## Remarks
  * For some C++ compilers the function expm1 is not available. A
    possible workaround is to use the expm1 from the boost libraries.
    Therefore, add the line #include <boost/math/special_functions/expm1.hpp\> 
    to the header file VineCopulaCPP_header.hpp. Then substitute expm1 by boost::math::expm1.
  * The function lgamma is also not available in some C++ compilers. A
    possible workaround is to use the lgamma function from the boost libraries.
    Therefore, add the line #include <boost/math/special_functions/gamma.hpp\> 
    to the header file VineCopulaCPP_header.hpp. Then substitute lgamma
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

# Parallel computing within the VineCopulaMatlab toolbox using OpenMP
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