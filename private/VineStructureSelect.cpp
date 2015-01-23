#include "VineCPP_header.hpp"
#include <matrix.h>
#include <mex.h>
    
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *familyset, *families, *rotations, *Thetas, *structure, *structuringrule;
    double *Type;
    int type, m, StructuringRule;
    unsigned int n;
    int i, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    d = (int)mxGetN(prhs[1]);
    
//associate inputs
    Type = mxGetPr(prhs[0]);
    type = (int) *Type;
    U = mxGetPr(prhs[1]);
    structuringrule = mxGetPr(prhs[2]);
    StructuringRule = (int) *structuringrule;
    m = (int)mxGetM(prhs[3]);
    familyset = mxGetPr(prhs[3]);
    
    std::vector<double> thetas;
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(1,d,mxREAL);
    structure = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,d*(d-1)/2,mxREAL);
    families = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,d*(d-1)/2,mxREAL);
    rotations = mxGetPr(plhs[2]);
    
    VineCopulaStructureSelect(type, structure, families, rotations, thetas, U, d, n, StructuringRule, familyset, m);
    
    int NumbParams = thetas.size();
    plhs[3] = mxCreateDoubleMatrix(1,NumbParams,mxREAL);
    Thetas = mxGetPr(plhs[3]);
    
    for (i=0;i<NumbParams;i++)
        Thetas[i] = thetas[i];
    
    return;
    
}
