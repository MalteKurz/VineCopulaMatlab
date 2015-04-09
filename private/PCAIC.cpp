#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *AIC, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    
//associate inputs    
    Family = mxGetPr(prhs[0]);
    family = (int) *Family;
    
// Depending on the family associate outputs
    switch(family){
        case 18:
        {
            plhs[1] = mxCreateDoubleMatrix(1,3,mxREAL);
            break;
        }
        case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
        {
            plhs[1] = mxCreateDoubleMatrix(1,2,mxREAL);
            break;
        }
        default:
        {
            plhs[1] = mxCreateDoubleScalar(mxREAL);
        }
    }
    
        
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    AIC = mxGetPr(plhs[0]);
    theta = mxGetPr(plhs[1]);
    
    if (nrhs==4 && !mxIsEmpty(prhs[3]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[1]);
        V_in = mxDuplicateArray(prhs[2]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[3]);
        rotation = (int) *Rotation;
        
        PairCopulaAIC(AIC,theta,family,rotation,U,V,n);
    }
    else
    {
        U = mxGetPr(prhs[1]);
        V = mxGetPr(prhs[2]);
        
        PairCopulaAIC(AIC,theta,family,U,V,n);
    }

    
    
    return;
    
}
