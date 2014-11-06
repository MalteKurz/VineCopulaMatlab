#include "VineCPP_header.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    double *U;
    double *N, *M;
    unsigned int i, n, m;
    

//figure out dimensions
    N = mxGetPr(prhs[0]);
    n = (unsigned int) *N;
    M = mxGetPr(prhs[1]);
    m = (unsigned int) *M;
    
//associate outputs
    plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL);
    U = mxGetPr(plhs[0]);
    
    boost::mt19937 gen;
    // Load the state
    std::ifstream fi(PathSeed);
    fi>>gen;
    fi.close();
    
    boost::uniform_01 <> URAND;
    
    boost::variate_generator <boost::mt19937&, boost::uniform_01 <> > RAND(gen,URAND);
    
    for (i=0;i<n*m;i++)
    {
        U[i] = RAND();
    }
    
    // Save the state
    std::ofstream fo(PathSeed,
            std::ios_base::out);
    fo<<gen;
    fo.close();
    
    return;
    
}
