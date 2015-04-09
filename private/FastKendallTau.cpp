#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *tau;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[0]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    
    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
    
    tau = mxGetPr(plhs[0]);
    
    *tau = SD_Kendall_Tau(U,V,n);

    return;
    
}
