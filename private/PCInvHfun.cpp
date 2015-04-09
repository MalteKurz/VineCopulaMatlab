#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *u, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n, N;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    N = (unsigned int)mxGetM(prhs[3]);
    
//associate inputs
    Family = mxGetPr(prhs[0]);
    family = (int) *Family;
    theta = mxGetPr(prhs[3]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    u = mxGetPr(plhs[0]);
    
    if (N == n)
    {
        if (nrhs==5 && !mxIsEmpty(prhs[4]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[1]);
            V_in = mxDuplicateArray(prhs[2]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[4]);
            rotation = (int) *Rotation;
            
            PairCopulaInvHfun_VecPar(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[1]);
            V = mxGetPr(prhs[2]);
            
            PairCopulaInvHfun_VecPar(family,theta,U,V,u,n);
        }
    }
    else
    {
        if (nrhs==5 && !mxIsEmpty(prhs[4]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[1]);
            V_in = mxDuplicateArray(prhs[2]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[4]);
            rotation = (int) *Rotation;
            
            PairCopulaInvHfun(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[1]);
            V = mxGetPr(prhs[2]);
            
            PairCopulaInvHfun(family,theta,U,V,u,n);
        }
    }
    
    
    return;
    
}
