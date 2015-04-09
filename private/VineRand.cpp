#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation;
    double *Type, *N, *D;
    int type;
    unsigned int n, d;
    
//associate inputs
    Type = mxGetPr(prhs[0]);
    type = (int) *Type;
    N = mxGetPr(prhs[1]);
    n = (unsigned int) *N;
    D = mxGetPr(prhs[2]);
    d = (unsigned int) *D;
    families = mxGetPr(prhs[3]);
    thetas = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,d,mxREAL);
    U = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        rotation = mxGetPr(prhs[5]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaRand(VinePtr,U,n);
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaRand(VinePtr,U,n);
    }



    
    return;
    
}
