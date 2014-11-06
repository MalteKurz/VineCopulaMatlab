# VineCopula (VineCopulaObject)

Constructor for objects of the VineCopula class

---
## Purpose
The "function" VineCopula can be used to specify
        members of the VineCopula class. At least the
        properties dimension, type and simplified have to be
        specified.


## Usage
        Three entries
               VineCopulaObject = VineCopula(dimension,type,simplified)
        Four entries
               VineCopulaObject = VineCopula(dimension,type,simplified,structure)
        Five entries
               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families)
        Six entries
               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families,parameters)
        Seven entries
               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families,parameters,condparameterfunctionals)

        Example 1: Three-dimensional non-simplified C-Vine
        copula (in three different but equivalent forms)
               VineCopulaObject = VineCopula(3,'C-Vine',[0],[1 2 3],{'Clayton','Gumbel','Frank'},{1.2, 3,[]},{@(u) (4.*u-2).^3})
               VineCopulaObject = VineCopula(3,'C-Vine',[0],[1 2 3],[7, 11, 9],{1.2, 3,[]},{@(u) (4.*u-2).^3})
               VineCopulaObject = VineCopula(3,'C-Vine',0,[1 2 3],{'Clayton','Gumbel','Frank'},[1.2, 3],{@(u) (4.*u-2).^3})
        Example 2: Three-dimensional simplified C-Vine copula
        (in three different but equivalent forms)
               VineCopulaObject = VineCopula(3,'C-Vine',[1],[1 2 3],{'Clayton','Gumbel','Indep'},{1.2, 3,[]})
               VineCopulaObject = VineCopula(3,'C-Vine',[1],[1 2 3],[7, 11, 0],{1.2, 3,[]})
               VineCopulaObject = VineCopula(3,'C-Vine',1,[1 2 3],{'Clayton','Gumbel','Indep'},[1.2, 3])

## References (for the examples):
Acar, E. F., C. Genest and J. Neslehová (2012), "Beyond
      simplified pair-copula constructions", Journal of
      Multivariate Analysis 110, pp. 74-90.
