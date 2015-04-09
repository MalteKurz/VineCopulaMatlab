#include "VineCopulaCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  boost::mt19937 gen;
  gen.seed(time(0));
  
  // Save the state
  std::ofstream fo(PathSeed, 
		   std::ios_base::out);
  fo<<gen;
  fo.close();
  
  return;
  
}
