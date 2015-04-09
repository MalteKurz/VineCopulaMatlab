#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation, *CLL;
    double *Type, *CutOffTree;
    int type, cutofftree;
    unsigned int n, d;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    d = (unsigned int)mxGetN(prhs[1]);
    
//associate inputs
    Type = mxGetPr(prhs[0]);
    type = (int) *Type;
    families = mxGetPr(prhs[2]);
    thetas = mxGetPr(prhs[3]);
    U = mxGetPr(prhs[1]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    CLL = mxGetPr(plhs[0]);
    
    //associate pointers and handle the rotations
    if (nrhs==6)
    {
        if (!mxIsEmpty(prhs[4]))
        {
        rotation = mxGetPr(prhs[4]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        
        CutOffTree = mxGetPr(prhs[5]);
        cutofftree = (int) *CutOffTree;
        *CLL = VineCopulaNegLL(VinePtr, U, cutofftree, n);
        }
        else
        {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        
        CutOffTree = mxGetPr(prhs[5]);
        cutofftree = (int) *CutOffTree;
        *CLL = VineCopulaNegLL(VinePtr, U, cutofftree, n);
        }
    }
    else
    {
        if (nrhs==5 && !mxIsEmpty(prhs[4]))
        {
            rotation = mxGetPr(prhs[4]);
            
            VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
            VineCopula* VinePtr = &Vine;
            *CLL = VineCopulaNegLL(VinePtr, U, d-1, n);
        }
        else
        {
            VineCopula Vine = VineCopula(type,d,families,thetas);
            VineCopula* VinePtr = &Vine;
            *CLL = VineCopulaNegLL(VinePtr, U, d-1, n);
        }
    }



    
    return;
    
}
