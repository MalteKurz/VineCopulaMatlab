#include "VineCPP_header.hpp"

bool PairCopulaIndepTest(double *U,double *V,unsigned int n, double *tau)
{
    *tau = SD_Kendall_Tau(U,V,n);
    
    double N1,N2;
    
    N1 = (double) 9*n*(n-1);
    N2 = (double) 2*(2*n+5);
    
    double TestStat = std::sqrt(N1/N2)*std::abs(*tau);
    
    boost::math::normal dist(0,1);

    double pValue = 2*(1-boost::math::cdf(dist, TestStat));
    
    //double pValue = 2*(1-gsl_cdf_ugaussian_P(TestStat));
    
    
    if (pValue>0.05)
            return true;
    else
        return false;
    
}

void PairCopulaIndepTest(double *pValue, double *U,double *V,unsigned int n)
{
    double tau = SD_Kendall_Tau(U,V,n);

    double N1,N2;
    
    N1 = (double) 9*n*(n-1);
    N2 = (double) 2*(2*n+5);
    
    double TestStat = std::sqrt(N1/N2)*std::abs(tau);
    
    boost::math::normal dist(0,1);

    *pValue = 2*(1-boost::math::cdf(dist, TestStat));
    //*pValue = 2*(1-gsl_cdf_ugaussian_P(TestStat));
    
}

void PairCopulaIndepTest(double *pValue, double *U,double *V,unsigned int n,double *TestStat)
{
    double tau = SD_Kendall_Tau(U,V,n);

    double N1,N2;
    
    N1 = (double) 9*n*(n-1);
    N2 = (double) 2*(2*n+5);
    
    *TestStat = std::sqrt(N1/N2)*std::abs(tau);
    
    boost::math::normal dist(0,1);

    *pValue = 2*(1-boost::math::cdf(dist, *TestStat));
    //*pValue = 2*(1-gsl_cdf_ugaussian_P(*TestStat));
    
}