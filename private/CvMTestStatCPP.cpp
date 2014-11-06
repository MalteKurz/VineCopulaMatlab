#include "VineCPP_header.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nlhs > 1 ) {
/*
 * This is the case, where only the value of the test statistic should be
 * given back to MATLAB.
 */
//declare variables
        const int *dims1, *dims2;
        double *U, *V, *TestStat, *V1, *V2, *V1dot, *V2dot, *V1dotdot, *V2dotdot;
        int d1, d2, n;
        int i,j,l;
        
//associate inputs
    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
        
//figure out dimensions
        dims1 = mxGetDimensions(prhs[0]);
        dims2 = mxGetDimensions(prhs[1]);
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
        const int *dims1, *dims2;
        double *U, *V, *TestStat;
        int d1, d2, n;
        int i,j,l;
        
//associate inputs
    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
        
//figure out dimensions
        dims1 = mxGetDimensions(prhs[0]);
        dims2 = mxGetDimensions(prhs[1]);
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
