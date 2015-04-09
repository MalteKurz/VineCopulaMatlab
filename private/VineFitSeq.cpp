#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation, *CLL;
    double *Type;
    int type, J=0;
    unsigned int i, n, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    d = (unsigned int)mxGetN(prhs[1]);
    
//associate inputs
    Type = mxGetPr(prhs[0]);
    type = (int) *Type;
    families = mxGetPr(prhs[2]);
    U = mxGetPr(prhs[1]);
    
    for (i=0;i<(d-1)*d/2;i++)
    {
        switch((int) families[i]){
            case 0:
            {
                break;
            }
            case 18:
            {
                J += 3;
                break;
            }
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                J += 2;
                break;
            }
            default:
            {
                J++;
            }
            
        }
    }
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(1,J,mxREAL);
    thetas = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleScalar(mxREAL);
    CLL = mxGetPr(plhs[1]);
    
    if (nrhs==4 && !mxIsEmpty(prhs[3]))
    {
        rotation = mxGetPr(prhs[3]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFitSeq(VinePtr, CLL, U, n);
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFitSeq(VinePtr, CLL, U, n);
    }



    
    return;
    
}
