#include "VineCPP_header.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *p, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    
//associate inputs
    Family = mxGetPr(prhs[0]);
    family = (int) *Family;
    theta = mxGetPr(prhs[3]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    p = mxGetPr(plhs[0]);
    
    if (nrhs==5 && !mxIsEmpty(prhs[4]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[1]);
        V_in = mxDuplicateArray(prhs[2]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[4]);
        rotation = (int) *Rotation;
        
        PairCopulaCDF(family,rotation,theta,U,V,p,n);
    }
    else
    {
        U = mxGetPr(prhs[1]);
        V = mxGetPr(prhs[2]);
        
        PairCopulaCDF(family,theta,U,V,p,n);
    }
    
    
    return;
    
}
