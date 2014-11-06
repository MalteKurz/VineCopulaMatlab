#Sim (VineCopulaObject)

Generating (pseudo-)random observations from a vine copula

---
## Purpose
The function draws N pseudo-random d-dimensional
        tuples from a vine copula. The function can be used to
        sample from arbitrarily dimensional C- and D-Vines of
        the simplified and non-simplified form. In the case of
        a non-simplified C- or D-Vine, the parameters of the
        conditional copulas have to be given as functionals of
        the variables on which the conditioning is done.


## Usage
       Simulate N observations from a vine copula
                      U = Sim(VineCopulaObject,N)


## Inputs
        VineCopulaObject= An object from the class VineCopula.
        N               = The number of observations that
                          should be simulated from the vine
                          copula.


## Outputs
        U               = A (N x d) Matrix of simulated tuples
                          from the specified vine copula, where
                          every row is a simulated d-dimensional
                          tuple from the d-dimensional vine
                          copula.
