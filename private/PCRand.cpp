#include "VineCPP_header.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *Rotation;
    double *Family, *N;
    int family, rotation;
    unsigned int i, n;
    
    
//associate inputs
    Family = mxGetPr(prhs[0]);
    family = (int) *Family;
    N = mxGetPr(prhs[1]);
    n = (unsigned int) *N;
    theta = mxGetPr(prhs[2]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    U = mxGetPr(plhs[0]);
    V = mxGetPr(plhs[1]);
    
    if (nrhs==4 && !mxIsEmpty(prhs[3]))
    {
        Rotation = mxGetPr(prhs[3]);
        rotation = (int) *Rotation;
        
        PairCopulaRand(family,rotation,theta,U,V,n);
    }
    else
    {
        PairCopulaRand(family,theta,U,V,n);
    }
    
    
    return;
    
}
