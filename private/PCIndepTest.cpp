#include "VineCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *pValue;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[0]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    pValue = mxGetPr(plhs[0]);
    
//associate inputs
    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
    
    if (nlhs==1)
    {
        PairCopulaIndepTest(pValue,U,V,n);
    }
    else
    {
        double *TestStat;
        plhs[1] = mxCreateDoubleScalar(mxREAL);
        TestStat = mxGetPr(plhs[1]);
        PairCopulaIndepTest(pValue,U,V,n,TestStat);
    }

    return;
    
}
