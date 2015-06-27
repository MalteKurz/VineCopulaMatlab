#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

// FunctionID 1: PCAIC
void PCAIC(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *AIC, *Rotation, *bounds;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[3]);
    
//associate inputs    
    bounds = mxGetPr(prhs[1]);
    Family = mxGetPr(prhs[2]);
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
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[3]);
        V_in = mxDuplicateArray(prhs[4]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaAIC(bounds,AIC,theta,family,rotation,U,V,n);
    }
    else
    {
        U = mxGetPr(prhs[3]);
        V = mxGetPr(prhs[4]);
        
        PairCopulaAIC(bounds,AIC,theta,family,U,V,n);
    }

    
    
    return;
    
}


// FunctionID 2: PCCDF
void PCCDF(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *p, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    p = mxGetPr(plhs[0]);
    
    if (nrhs==7 && !mxIsEmpty(prhs[4]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[2]);
        V_in = mxDuplicateArray(prhs[3]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[4]);
        rotation = (int) *Rotation;
        
        PairCopulaCDF(family,rotation,theta,U,V,p,n);
    }
    else
    {
        U = mxGetPr(prhs[2]);
        V = mxGetPr(prhs[3]);
        
        PairCopulaCDF(family,theta,U,V,p,n);
    }
    
    
    return;
    
}


// FunctionID 3: PCFit
void PCFit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *Rotation, *bounds;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[3]);
    
//associate inputs
    bounds = mxGetPr(prhs[1]);
    Family = mxGetPr(prhs[2]);
    family = (int) *Family;
    
// Depending on the family associate outputs and upper and lower bounds
    switch(family){
        case 18:
        {
            plhs[0] = mxCreateDoubleMatrix(1,3,mxREAL);
            break;
        }
        case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
        {
            plhs[0] = mxCreateDoubleMatrix(1,2,mxREAL);
            break;
        }
        default:
        {
            plhs[0] = mxCreateDoubleScalar(mxREAL);
        }
    }
    
//associate outputs
    theta = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[3]);
        V_in = mxDuplicateArray(prhs[4]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaFit(bounds,theta,family,rotation,U,V,n);
    }
    else
    {
        U = mxGetPr(prhs[3]);
        V = mxGetPr(prhs[4]);
        
        PairCopulaFit(bounds,theta,family,U,V,n);
    }
    
    

    
    return;
    
}

// FunctionID 4: PCHfun
void PCHfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *u, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    u = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[2]);
        V_in = mxDuplicateArray(prhs[3]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaHfun(family,rotation,theta,U,V,u,n);
        
    }
    else
    {
        U = mxGetPr(prhs[2]);
        V = mxGetPr(prhs[3]);
        
        PairCopulaHfun(family,theta,U,V,u,n);
    }
    
    
    return;
    
}

// FunctionID 5: PCIndepTest
void PCIndepTest(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *pValue;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    pValue = mxGetPr(plhs[0]);
    
//associate inputs
    U = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
    
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

// FunctionID 6: PCInvHfun
void PCInvHfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *u, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n, N;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    N = (unsigned int)mxGetM(prhs[4]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    u = mxGetPr(plhs[0]);
    
    if (N == n)
    {
        if (nrhs==6 && !mxIsEmpty(prhs[5]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[2]);
            V_in = mxDuplicateArray(prhs[3]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[5]);
            rotation = (int) *Rotation;
            
            PairCopulaInvHfun_VecPar(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[2]);
            V = mxGetPr(prhs[3]);
            
            PairCopulaInvHfun_VecPar(family,theta,U,V,u,n);
        }
    }
    else
    {
        if (nrhs==6 && !mxIsEmpty(prhs[5]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[2]);
            V_in = mxDuplicateArray(prhs[3]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[5]);
            rotation = (int) *Rotation;
            
            PairCopulaInvHfun(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[2]);
            V = mxGetPr(prhs[3]);
            
            PairCopulaInvHfun(family,theta,U,V,u,n);
        }
    }
    
    
    return;
    
}

// FunctionID 7: PCInvVfun
void PCInvVfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *u, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n, N;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    N = (unsigned int)mxGetM(prhs[4]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    u = mxGetPr(plhs[0]);
    
    if (N == n)
    {
        if (nrhs==6 && !mxIsEmpty(prhs[5]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[2]);
            V_in = mxDuplicateArray(prhs[3]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[5]);
            rotation = (int) *Rotation;
            
            PairCopulaInvVfun_VecPar(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[2]);
            V = mxGetPr(prhs[3]);
            
            PairCopulaInvVfun_VecPar(family,theta,U,V,u,n);
        }
    }
    else
    {
        if (nrhs==6 && !mxIsEmpty(prhs[5]))
        {
            mxArray *U_in, *V_in;
            U_in = mxDuplicateArray(prhs[2]);
            V_in = mxDuplicateArray(prhs[3]);
            U = mxGetPr(U_in);
            V = mxGetPr(V_in);
            Rotation = mxGetPr(prhs[5]);
            rotation = (int) *Rotation;
            
            PairCopulaInvVfun(family,rotation,theta,U,V,u,n);
            
        }
        else
        {
            U = mxGetPr(prhs[2]);
            V = mxGetPr(prhs[3]);
            
            PairCopulaInvVfun(family,theta,U,V,u,n);
        }
    }
    
    
    return;
    
}

// FunctionID 8: PCNegLL
void PCNegLL(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *CLL, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    CLL = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[2]);
        V_in = mxDuplicateArray(prhs[3]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        Rotate_Obs(U,V,rotation,n);
    }
    else
    {
        U = mxGetPr(prhs[2]);
        V = mxGetPr(prhs[3]);
    }
    
    *CLL = PairCopulaNegLL(family,theta,U,V,n);
    
    return;
    
}

// FunctionID 9: PCPDF
void PCPDF(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *p, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    p = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[2]);
        V_in = mxDuplicateArray(prhs[3]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaPDF(family,rotation,theta,U,V,p,n);
    }
    else
    {
        U = mxGetPr(prhs[2]);
        V = mxGetPr(prhs[3]);
        
        PairCopulaPDF(family,theta,U,V,p,n);
    }
    
    
    return;
    
}

// FunctionID 10: PCRand
void PCRand(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *Rotation;
    double *Family, *N;
    int family, rotation;
    unsigned int n;
    
    
//associate inputs
    Family = mxGetPr(prhs[2]);
    family = (int) *Family;
    N = mxGetPr(prhs[3]);
    n = (unsigned int) *N;
    theta = mxGetPr(prhs[4]);
    
    
// Load the state of the seed
    double *OldState, *StateExport;
    int i;
    OldState = mxGetPr(prhs[1]);
    
    std::vector<unsigned int> StateImport(624);
    for (i=0;i<624;i++)
    {
        StateImport[i] = (unsigned int) OldState[i];
    }
    
//associate outputs
    plhs[1] = mxCreateDoubleMatrix(n,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
    U = mxGetPr(plhs[1]);
    V = mxGetPr(plhs[2]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaRand(StateImport,family,rotation,theta,U,V,n);
    }
    else
    {
        PairCopulaRand(StateImport,family,theta,U,V,n);
    }
    
    plhs[0] = mxCreateDoubleMatrix(624, 1, mxREAL);
    StateExport = mxGetPr(plhs[0]);
    for (i=0;i<624;i++)
    {
        StateExport[i] = (double) StateImport[i];
    }
    
    
    return;
    
}

// FunctionID 11: PCSelect
void PCSelect(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *Theta, *Rotation, *Family,*familyset, *bounds;
    int family, rotation;
    unsigned int n;
    int m;
    std::vector<double> theta(3);
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    m = (int)mxGetM(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    Family = mxGetPr(plhs[0]);
    
//associate inputs
    bounds = mxGetPr(prhs[1]);
    mxArray *U_in, *V_in;
    U_in = mxDuplicateArray(prhs[2]);
    V_in = mxDuplicateArray(prhs[3]);
    U = mxGetPr(U_in);
    V = mxGetPr(V_in);
    familyset = mxGetPr(prhs[4]);
    
    PairCopulaSelect(bounds,&family,&theta[0],&rotation,U,V,n,familyset,m);
    
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

// FunctionID 12: PCVfun
void PCVfun(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *theta, *u, *Rotation;
    double *Family;
    int family, rotation;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    
//associate inputs
    Family = mxGetPr(prhs[1]);
    family = (int) *Family;
    theta = mxGetPr(prhs[4]);
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
    u = mxGetPr(plhs[0]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        mxArray *U_in, *V_in;
        U_in = mxDuplicateArray(prhs[2]);
        V_in = mxDuplicateArray(prhs[3]);
        U = mxGetPr(U_in);
        V = mxGetPr(V_in);
        Rotation = mxGetPr(prhs[5]);
        rotation = (int) *Rotation;
        
        PairCopulaVfun(family,rotation,theta,U,V,u,n);
    }
    else
    {
        U = mxGetPr(prhs[2]);
        V = mxGetPr(prhs[3]);
        
        PairCopulaVfun(family,theta,U,V,u,n);
    }
    
    
    return;
    
}

// FunctionID 101: VineFit
void VineFit(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *thetas0, *rotation, *CLL, *bounds;
    double *Type;
    int type, J=0;
    unsigned int i, n, d;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[3]);
    d = (unsigned int)mxGetN(prhs[3]);
    
//associate inputs
    bounds = mxGetPr(prhs[1]);
    Type = mxGetPr(prhs[2]);
    type = (int) *Type;
    families = mxGetPr(prhs[4]);
    U = mxGetPr(prhs[3]);
    
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
    plhs[1] = mxCreateDoubleMatrix(1,2,mxREAL);
    CLL = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,J,mxREAL);
    thetas0 = mxGetPr(plhs[2]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        rotation = mxGetPr(prhs[5]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFit(bounds, VinePtr, CLL, thetas0, U, n);
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFit(bounds, VinePtr, CLL, thetas0, U, n);
    }



    
    return;
    
}

// FunctionID 102: VineFitSeq
void VineFitSeq(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation, *CLL, *bounds;
    double *Type;
    int type, J=0;
    unsigned int i, n, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[3]);
    d = (unsigned int)mxGetN(prhs[3]);
    
//associate inputs
    bounds = mxGetPr(prhs[1]);
    Type = mxGetPr(prhs[2]);
    type = (int) *Type;
    families = mxGetPr(prhs[4]);
    U = mxGetPr(prhs[3]);
    
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
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        rotation = mxGetPr(prhs[5]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFitSeq(bounds,VinePtr, CLL, U, n);
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaFitSeq(bounds,VinePtr, CLL, U, n);
    }



    
    return;
    
}

// FunctionID 103: VineGetPseudoObs
void VineGetPseudoObs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *V1, *V2, *families, *thetas, *rotation;
    double *Type;
    int type;
    unsigned int n, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    d = (unsigned int)mxGetN(prhs[2]);
    
//associate inputs
    Type = mxGetPr(prhs[1]);
    type = (int) *Type;
    families = mxGetPr(prhs[3]);
    thetas = mxGetPr(prhs[4]);
    U = mxGetPr(prhs[2]);
    
    if (nrhs==6 && !mxIsEmpty(prhs[5]))
    {
        rotation = mxGetPr(prhs[5]);
        
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

// FunctionID 104: VineNegLL
void VineNegLL(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation, *CLL;
    double *Type, *CutOffTree;
    int type, cutofftree;
    unsigned int n, d;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[2]);
    d = (unsigned int)mxGetN(prhs[2]);
    
//associate inputs
    Type = mxGetPr(prhs[1]);
    type = (int) *Type;
    families = mxGetPr(prhs[3]);
    thetas = mxGetPr(prhs[4]);
    U = mxGetPr(prhs[2]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    CLL = mxGetPr(plhs[0]);
    
    //associate pointers and handle the rotations
    if (nrhs==7)
    {
        if (!mxIsEmpty(prhs[5]))
        {
        rotation = mxGetPr(prhs[5]);
        
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        
        CutOffTree = mxGetPr(prhs[6]);
        cutofftree = (int) *CutOffTree;
        *CLL = VineCopulaNegLL(VinePtr, U, cutofftree, n);
        }
        else
        {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        
        CutOffTree = mxGetPr(prhs[6]);
        cutofftree = (int) *CutOffTree;
        *CLL = VineCopulaNegLL(VinePtr, U, cutofftree, n);
        }
    }
    else
    {
        if (nrhs==6 && !mxIsEmpty(prhs[5]))
        {
            rotation = mxGetPr(prhs[5]);
            
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

// FunctionID 105: VineRand
void VineRand(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *families, *thetas, *rotation;
    double *Type, *N, *D;
    int type;
    unsigned int n, d;
    
//associate inputs
    Type = mxGetPr(prhs[2]);
    type = (int) *Type;
    N = mxGetPr(prhs[3]);
    n = (unsigned int) *N;
    D = mxGetPr(prhs[4]);
    d = (unsigned int) *D;
    families = mxGetPr(prhs[5]);
    thetas = mxGetPr(prhs[6]);
    
    
// Load the state of the seed
    double *OldState, *StateExport;
    int i;
    OldState = mxGetPr(prhs[1]);
    
    std::vector<unsigned int> StateImport(624);
    for (i=0;i<624;i++)
    {
        StateImport[i] = (unsigned int) OldState[i];
    }
    
//associate outputs
    plhs[1] = mxCreateDoubleMatrix(n,d,mxREAL);
    U = mxGetPr(plhs[1]);
    
    if (nrhs==8 && !mxIsEmpty(prhs[7]))
    {
        rotation = mxGetPr(prhs[7]);
         
        VineCopula Vine = VineCopula(type,d,families,rotation,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaRand(StateImport,VinePtr,U,n);
    }
    else
    {
        VineCopula Vine = VineCopula(type,d,families,thetas);
        VineCopula* VinePtr = &Vine;
        VineCopulaRand(StateImport,VinePtr,U,n);
    }
    
    plhs[0] = mxCreateDoubleMatrix(624, 1, mxREAL);
    StateExport = mxGetPr(plhs[0]);
    for (i=0;i<624;i++)
    {
        StateExport[i] = (double) StateImport[i];
    }
    
    
    return;
    
}

// FunctionID 106: VineStructureSelect
void VineStructureSelect(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *familyset, *families, *rotations, *Thetas, *structure, *structuringrule, *bounds;
    double *Type;
    int type, m, StructuringRule;
    unsigned int n;
    int i, d;
    

//figure out dimensions
    n = (unsigned int)mxGetM(prhs[3]);
    d = (int)mxGetN(prhs[3]);
    
//associate inputs
    bounds = mxGetPr(prhs[1]);
    Type = mxGetPr(prhs[2]);
    type = (int) *Type;
    U = mxGetPr(prhs[3]);
    structuringrule = mxGetPr(prhs[4]);
    StructuringRule = (int) *structuringrule;
    m = (int)mxGetM(prhs[5]);
    familyset = mxGetPr(prhs[5]);
    
    std::vector<double> thetas;
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(1,d,mxREAL);
    structure = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1,d*(d-1)/2,mxREAL);
    families = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,d*(d-1)/2,mxREAL);
    rotations = mxGetPr(plhs[2]);
    
    VineCopulaStructureSelect(bounds,type, structure, families, rotations, thetas, U, d, n, StructuringRule, familyset, m);
    
    int NumbParams = thetas.size();
    plhs[3] = mxCreateDoubleMatrix(1,NumbParams,mxREAL);
    Thetas = mxGetPr(plhs[3]);
    
    for (i=0;i<NumbParams;i++)
        Thetas[i] = thetas[i];
    
    return;
    
}

// FunctionID 1001: CvMTestStatCPP
void CvMTestStatCPP(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs > 1 ) {
/*
 * This is the case, where only the value of the test statistic should be
 * given back to MATLAB.
 */

//declare variables
        const mwSize *dims1, *dims2;
        double *U, *V, *TestStat, *V1, *V2, *V1dot, *V2dot, *V1dotdot, *V2dotdot;
        int d1, d2, n;
        int i,j,l;
        
//associate inputs
    U = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
        
//figure out dimensions
        dims1 = mxGetDimensions(prhs[1]);
        dims2 = mxGetDimensions(prhs[2]);
        n = (int)dims1[0];
        d1 = (int)dims1[1];
        d2 = (int)dims2[1];
        
//associate outputs
        plhs[0] = mxCreateDoubleScalar(mxREAL);
        plhs[1] = mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[2] = mxCreateDoubleMatrix(n,n,mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1,n,mxREAL);
        plhs[4] = mxCreateDoubleMatrix(1,n,mxREAL);
        plhs[5] = mxCreateDoubleScalar(mxREAL);
        plhs[6] = mxCreateDoubleScalar(mxREAL);
        
        TestStat = mxGetPr(plhs[0]);
        V1 = mxGetPr(plhs[1]);
        V2 = mxGetPr(plhs[2]);
        V1dot = mxGetPr(plhs[3]);
        V2dot = mxGetPr(plhs[4]);
        V1dotdot = mxGetPr(plhs[5]);
        V2dotdot = mxGetPr(plhs[6]);
        
//doing the computations
        *V1dotdot = 0;
        *V2dotdot = 0;
        *TestStat = 0;
        
        for (i=0;i<n;i++)
        {
            *(V1+i*n+i) = 1-U[i];
            if (d1>1)
            {
                for (l=1;l<d1;l++)
                {
                    *(V1+i*n+i) *= (1-U[l*n+i]);
                }
            }
            for (j=i+1;j<n;j++)
            {
                *(V1+j*n+i) = 1-max(U[i],U[j]);
                *(V1+i*n+j) = *(V1+j*n+i);
            }
            if (d1>1)
            {
                for (l=1;l<d1;l++)
                {
                    for (j=i+1;j<n;j++)
                    {
                        *(V1+j*n+i) *= (1-max(U[l*n+i],U[l*n+j]));
                        *(V1+i*n+j) = *(V1+j*n+i);
                    }
                }
            }
        }
        
        for (i=0;i<n;i++)
        {
            *(V2+i*n+i) = 1-V[i];
            if (d2>1)
            {
                for (l=1;l<d2;l++)
                {
                    *(V2+i*n+i) *= (1-V[l*n+i]);
                }
            }
            for (j=i+1;j<n;j++)
            {
                *(V2+j*n+i) = 1-max(V[i],V[j]);
                *(V2+i*n+j) = *(V2+j*n+i);
            }
            if (d2>1)
            {
                for (l=1;l<d2;l++)
                {
                    for (j=i+1;j<n;j++)
                    {
                        *(V2+j*n+i) *= (1-max(V[l*n+i],V[l*n+j]));
                        *(V2+i*n+j) = *(V2+j*n+i);
                    }
                }
            }
        }
        
        for (j=0;j<n;j++)
        {
            for (i=0;i<n;i++)
            {
                *(V1dot+j) += (*(V1+j*n+i)/(double)n);
                *(V2dot+j) += (*(V2+j*n+i)/(double)n);
                *TestStat += ((*(V1+j*n+i) * *(V2+j*n+i))/(double)n);
            }
            *V1dotdot += (*(V1dot+j)/(double)n);
            *V2dotdot += (*(V2dot+j)/(double)n);
            *TestStat -= (2*(*(V1dot+j) * *(V2dot+j)));
        }
        
        *TestStat += ((double)n * *V1dotdot * *V2dotdot);
    }
    else {
 /*
 * This is the case, where the value of the test statistic should be given
 * back to MATLAB, but also the values of the matrices V1 and V2, the
 * vectors V1dot and V2dot and the scalar values V1dotdot and V2dotdot, 
 * which altogether can be used for computing the boostrapped values of the
 * test statistics in the case of using the multiplier bootstrap approach.
 */
        //declare variables
        const mwSize *dims1, *dims2;
        double *U, *V, *TestStat;
        int d1, d2, n;
        int i,j,l;
        
//associate inputs
    U = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
        
//figure out dimensions
        dims1 = mxGetDimensions(prhs[1]);
        dims2 = mxGetDimensions(prhs[2]);
        n = (int)dims1[0];
        d1 = (int)dims1[1];
        d2 = (int)dims2[1];
        
//associate outputs
        plhs[0] = mxCreateDoubleScalar(mxREAL);
        TestStat = mxGetPr(plhs[0]);
        
//doing the computations
        std::vector<double> V1(n*(n+1)/2);
        std::vector<double> V2(n*(n+1)/2);
        std::vector<double> V1dot(n);
        std::vector<double> V2dot(n);
        double V1dotdot = 0;
        double V2dotdot = 0;
        *TestStat = 0;
        
        for (i=0;i<n;i++)
        {
            V1[IndTria(i,i,n)] = 1-U[i];
            if (d1>1)
            {
                for (l=1;l<d1;l++)
                {
                    V1[IndTria(i,i,n)] *= (1-U[l*n+i]);
                }
            }
            for (j=i+1;j<n;j++)
            {
                V1[IndTria(i,j,n)] = 1-max(U[i],U[j]);
            }
            if (d1>1)
            {
                for (l=1;l<d1;l++)
                {
                    for (j=i+1;j<n;j++)
                    {
                        V1[IndTria(i,j,n)] *= (1-max(U[l*n+i],U[l*n+j]));
                    }
                }
            }
        }
        
        for (i=0;i<n;i++)
        {
            V2[IndTria(i,i,n)] = 1-V[i];
            if (d2>1)
            {
                for (l=1;l<d2;l++)
                {
                    V2[IndTria(i,i,n)] *= (1-V[l*n+i]);
                }
            }
            for (j=i+1;j<n;j++)
            {
                V2[IndTria(i,j,n)] = 1-max(V[i],V[j]);
            }
            if (d2>1)
            {
                for (l=1;l<d2;l++)
                {
                    for (j=i+1;j<n;j++)
                    {
                        V2[IndTria(i,j,n)] *= (1-max(V[l*n+i],V[l*n+j]));
                    }
                }
            }
        }
        
        for (j=0;j<n;j++)
        {
            for (i=0;i<n;i++)
            {
                V1dot[j] += (V1[IndTria(i,j,n)]/(double)n);
                V2dot[j] += (V2[IndTria(i,j,n)]/(double)n);
                *TestStat += ((V1[IndTria(i,j,n)]*V2[IndTria(i,j,n)])/(double)n);
            }
            V1dotdot += (V1dot[j]/(double)n);
            V2dotdot += (V2dot[j]/(double)n);
            *TestStat -= (2*(V1dot[j]*V2dot[j]));
        }
        
        *TestStat += ((double)n*V1dotdot*V2dotdot);
    }
    
    return;
}

// FunctionID 1002: FastKendallTau
void FastKendallTau(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U, *V, *tau;
    unsigned int n;
    
//figure out dimensions
    n = (unsigned int)mxGetM(prhs[1]);
    
//associate outputs
    plhs[0] = mxCreateDoubleScalar(mxREAL);
    
    U = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
    
    tau = mxGetPr(plhs[0]);
    
    *tau = SD_Kendall_Tau(U,V,n);

    return;
    
}


// FunctionID 1003: RandNormal
void RandNormal(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U;
    double *N, *M, *OldState, *StateExport;
    unsigned int i, n, m;
        
// Load the state of the seed
    boost::mt19937 gen;
    OldState = mxGetPr(prhs[1]);
    
    std::vector<unsigned int> StateImport(624);
    for (i=0;i<624;i++)
    {
        StateImport[i] = (unsigned int) OldState[i];
    }
    
    std::stringstream RNG_State_In;
    std::copy(StateImport.begin(), StateImport.end(), std::ostream_iterator<unsigned int>(RNG_State_In, " "));
    RNG_State_In>>gen;

//figure out dimensions
    N = mxGetPr(prhs[2]);
    n = (unsigned int) *N;
    M = mxGetPr(prhs[3]);
    m = (unsigned int) *M;
    
//associate outputs
    plhs[1] = mxCreateDoubleMatrix(n,m,mxREAL);
    U = mxGetPr(plhs[1]);
    
    
    boost::normal_distribution <> NRAND;
    
    boost::variate_generator<boost::mt19937& ,boost::normal_distribution <> > RAND(gen,NRAND);
    
    for (i=0;i<n*m;i++)
    {
        U[i] = RAND();
    }
    
    // Save the state
    double state;
    plhs[0] = mxCreateDoubleMatrix(624, 1, mxREAL);
    StateExport = mxGetPr(plhs[0]);
    std::stringstream RNG_State_Out;
    
    RNG_State_Out << gen;
    for (i=0;i<624;i++)
    {
        RNG_State_Out >> state;
        StateExport[i] = (double) state;
    }
    
    
    
    return;
    
}

// FunctionID 1004: RandUniform
void RandUniform(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U;
    double *N, *M, *OldState, *StateExport;
    unsigned int n, m;
    int i;
    
// Load the state of the seed
    boost::mt19937 gen;
    OldState = mxGetPr(prhs[1]);
    
    std::vector<unsigned int> StateImport(624);
    for (i=0;i<624;i++)
    {
        StateImport[i] = (unsigned int) OldState[i];
    }
    
    std::stringstream RNG_State_In;
    std::copy(StateImport.begin(), StateImport.end(), std::ostream_iterator<unsigned int>(RNG_State_In, " "));
    RNG_State_In>>gen;
    
//figure out dimensions
    N = mxGetPr(prhs[2]);
    n = (unsigned int) *N;
    M = mxGetPr(prhs[3]);
    m = (unsigned int) *M;
    
//associate outputs
    plhs[1] = mxCreateDoubleMatrix(n,m,mxREAL);
    U = mxGetPr(plhs[1]);
    
    
    boost::uniform_01 <> URAND;
    
    boost::variate_generator <boost::mt19937&, boost::uniform_01 <> > RAND(gen,URAND);
    
    for (i=0;i<n*m;i++)
    {
        U[i] = RAND();
    }
    
    // Save the state
    double state;
    plhs[0] = mxCreateDoubleMatrix(624, 1, mxREAL);
    StateExport = mxGetPr(plhs[0]);
    std::stringstream RNG_State_Out;
    
    RNG_State_Out << gen;
    for (i=0;i<624;i++)
    {
        RNG_State_Out >> state;
        StateExport[i] = (double) state;
    }
    
    return;
    
}

// FunctionID 1005: VineCopulaMatlabSetSeed
void VineCopulaMatlabSetSeed(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    boost::mt19937 gen;
    
    // Switch depending on the number of inputs
    if (nrhs==2 && !mxIsEmpty(prhs[1]))
    {    double *SeedInput;
         int IntSeed;
         
         SeedInput = mxGetPr(prhs[1]);
         IntSeed = (int) SeedInput[0];
         
         gen.seed(IntSeed);
    }
    else
    {
        gen.seed(time(0));
    }
    
//declare variables
    double *StateExport;
    unsigned int i;
    
//associate outputs
    std::stringstream RNG_State_Out;
    RNG_State_Out << gen;
    
    plhs[0] = mxCreateDoubleMatrix(624, 1, mxREAL);
    StateExport = mxGetPr(plhs[0]);
    
    double state;
    
    for (i=0;i<624;i++)
    {
        RNG_State_Out >> state;
        StateExport[i] = state;
    }

  return;
  
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *FunctionID;
    int functionID;
//associate inputs    
    FunctionID = mxGetPr(prhs[0]);
    functionID = (int) *FunctionID;
    
// Depending on the functionID call the corresponding function
    switch(functionID){
        case 1: // PCAIC
        {
            PCAIC(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 2: // PCCDF
        {
            PCCDF(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 3: // PCFit
        {
            PCFit(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 4: // PCHfun
        {
            PCHfun(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 5: // PCIndepTest
        {
            PCIndepTest(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 6: // PCInvHfun
        {
            PCInvHfun(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 7: // PCInvVfun
        {
            PCInvVfun(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 8: // PCNegLL
        {
            PCNegLL(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 9: // PCPDF
        {
            PCPDF(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 10: // PCRand
        {
            PCRand(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 11: // PCSelect
        {
            PCSelect(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 12: // PCVfun
        {
            PCVfun(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 101: // VineFit
        {
            VineFit(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 102: // VineFitSeq
        {
            VineFitSeq(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 103: // VineGetPseudoObs
        {
            VineGetPseudoObs(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 104: // VineNegLL
        {
            VineNegLL(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 105: // VineRand
        {
            VineRand(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 106: // VineStructureSelect
        {
            VineStructureSelect(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 1001: // CvMTestStatCPP
        {
            CvMTestStatCPP(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 1002: // FastKendallTau
        {
            FastKendallTau(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 1003: // RandNormal
        {
            RandNormal(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 1004: // RandUniform
        {
            RandUniform(nlhs, plhs, nrhs, prhs);
            break;
        }
        case 1005: // VineCopulaMatlabSetSeed
        {
            VineCopulaMatlabSetSeed(nlhs, plhs, nrhs, prhs);
            break;
        }
        default:
        {
            // Place an error here
        }
    }
    
    return;
    
}


