# SeqTestOnSimplified (VineCopulaObject)

Sequentially testing the simplified assumption for vine copulas

---
## Purpose
The function can be used to sequentially test the
        simplifying assumption for vine copulas. The procedure
        heavily relies on a stochastic representation of the
        simplifying assumption stated in Spanhel and Kurz
        (2014). The stochastic representation allows to test
        the simplifying assumption by using tests on vectorial
        independencies. Note, that the performed tests are 
        always based on the assumption that one is able to
        observe (pseudo-)observations from the conditional
        copulas that should be tested. Therefore, for
        interpreting the test results for a specific tree of
        the vine copula one needs to assume that the lower
        trees (including the assumption of unconditional
        copulas) are correctly specified.


## Usage
        Testing the simplifying assumption sequentially
         [pVals,TestStats,BootTestStats] = SeqTestOnSimplified(VineCopulaObject,data,N)
        Testing the simplifying assumption sequentially as
        goodness-of-fit test (i.e., without reestimating the
        whole vine copula for each tree)
         [pVals,TestStats,BootTestStats] = SeqTestOnSimplified(VineCopulaObject,data,N,true)


## Inputs
        VineCopulaObject      = An object from the class
                                VineCopula.
        data                  = A (n x d) dimensional vector
                                of values lying in [0,1] (the
                                observations).
        N                     = The number of boostrap
                                replications for the
                                multiplier bootstrap in the
                                vectorial independence test.
        GoF                   = A logical, which is by default
                                false. Then, before the tests
                                on vectorial independence are
                                applied in each tree, the
                                whole vine copula model is
                                reestimated up to this tree
                                and then the tests are
                                performed. Otherwise, i.e., if
                                it is applied as a goodness-
                                of-fit test, the joint
                                estimates given as inputs
                                through the VineCopula object
                                are used to obtain the
                                (pseudo-)observations in each
                                tree (without re-estimation).


## Outputs
        pVal                  = A vector of p-values of the
                                vectorial independence tests. 
                                Every entry corresponds to one
                                specific copula being part
                                of the whole vine copula.
        TestStat              = A vector of realized values
                                for the test statistics.
                                of the whole vine copula.
        BootTestStats         = A matrix of realized values
                                for the bootstrapped test
                                statistics.
