#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *Theta, *Rotation, *Family,*familyset;
    int family, rotation;
    unsigned int n;
    int m;
    std::vector<double> theta(3);
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[0]);
    m = (int)mxGetM(prhs[2]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    Family = mxGetPr(plhs[0]);
    
//associate inputs
    mxArray *U_in, *V_in;
    U_in = mxDuplicateArray(prhs[0]);
    V_in = mxDuplicateArray(prhs[1]);
    U = mxGetPr(U_in);
    V = mxGetPr(V_in);
    familyset = mxGetPr(prhs[2]);
    
    PairCopulaSelect(&family,&theta[0],&rotation,U,V,n,familyset,m);
    
    *Family = (double) family;
// Depending on the family associate outputs and upper and lower bounds
    switch(family){
        case 0:
        {
            plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
            plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
            break;
        }
        case 18:
        {
            plhs[1] = mxCreateDoubleMatrix(1,3,mxREAL);
            Theta = mxGetPr(plhs[1]);
            Theta[0] = theta[0];
            Theta[1] = theta[1];
            Theta[2] = theta[2];
            if (rotation==0)
            {
                plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
            }
            else
            {
                plhs[2] = mxCreateDoubleScalar(mxREAL);
                Rotation = mxGetPr(plhs[2]);
                *Rotation = (double) rotation;
            }
            break;
        }
        case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
        {
            plhs[1] = mxCreateDoubleMatrix(1,2,mxREAL);
            Theta = mxGetPr(plhs[1]);
            Theta[0] = theta[0];
            Theta[1] = theta[1];
            if (rotation==0)
            {
                plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
            }
            else
            {
                plhs[2] = mxCreateDoubleScalar(mxREAL);
                Rotation = mxGetPr(plhs[2]);
                *Rotation = (double) rotation;
            }
            break;
        }
        default:
        {
            plhs[1] = mxCreateDoubleScalar(mxREAL);
            Theta = mxGetPr(plhs[1]);
            Theta[0] = theta[0];
            if (rotation==0)
            {
                plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
            }
            else
            {
                plhs[2] = mxCreateDoubleScalar(mxREAL);
                Rotation = mxGetPr(plhs[2]);
                *Rotation = (double) rotation;
            }
        }
    }
    
    return;
    
}
