#include "VineCPP_header.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *V1, *V2, *families, *thetas, *rotation;
    double *Type;
    int type;
    unsigned int i, n, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    d = (unsigned int)mxGetN(prhs[1]);
    
//associate inputs
    Type = mxGetPr(prhs[0]);
    type = (int) *Type;
    families = mxGetPr(prhs[2]);
    thetas = mxGetPr(prhs[3]);
    U = mxGetPr(prhs[1]);
    
    if (nrhs==5 && !mxIsEmpty(prhs[4]))
    {
        rotation = mxGetPr(prhs[4]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        
        switch (type)
        {
            case 0:
            {
                plhs[0] = mxCreateDoubleMatrix(n,d*(d-1)/2-1,mxREAL);
                V = mxGetPr(plhs[0]);
                VineCopulaGetPseudoObs(VinePtr, U, V, n);
                break;
            }
            case 1:
            {
                plhs[0] = mxCreateDoubleMatrix(n,(d-1)*(d-2)/2,mxREAL);
                V1 = mxGetPr(plhs[0]);
                plhs[1] = mxCreateDoubleMatrix(n,(d-1)*(d-2)/2,mxREAL);
                V2 = mxGetPr(plhs[1]);
                VineCopulaGetPseudoObs(VinePtr, U, V1, V2, n);
            }
        }
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        
        switch (type)
        {
            case 0:
            {
                plhs[0] = mxCreateDoubleMatrix(n,d*(d-1)/2-1,mxREAL);
                V = mxGetPr(plhs[0]);
                VineCopulaGetPseudoObs(VinePtr, U, V, n);
                break;
            }
            case 1:
            {
                plhs[0] = mxCreateDoubleMatrix(n,(d-1)*(d-2)/2,mxREAL);
                V1 = mxGetPr(plhs[0]);
                plhs[1] = mxCreateDoubleMatrix(n,(d-1)*(d-2)/2,mxREAL);
                V2 = mxGetPr(plhs[1]);
                VineCopulaGetPseudoObs(VinePtr, U, V1, V2, n);
            }
        }
    }



    
    return;
    
}
