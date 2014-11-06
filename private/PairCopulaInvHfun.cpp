#include "VineCPP_header.hpp"

void PairCopulaInvHfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{  
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvHfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 180:
        {
            PairCopulaInvHfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            break;
        }
    }
    
    return;
    
}

void PairCopulaInvVfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvVfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvHfun(family, theta, V, U, u, n);
            break;
        }
        case 180:
        {
            PairCopulaInvVfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvHfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
    }
    
    return;
    
}

void PairCopulaInvHfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvHfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 180:
        {
            PairCopulaInvHfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            break;
        }
    }
    
    return;
    
}

void PairCopulaInvVfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvVfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvHfun(family, theta, V, U, u, n);
            break;
        }
        case 180:
        {
            PairCopulaInvVfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvHfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
    }
    
    return;
    
}

void PairCopulaInvHfun_VecPar(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvHfun_VecPar(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 180:
        {
            PairCopulaInvHfun_VecPar(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvVfun(family, theta, V, U, u, n);
            break;
        }
    }
    
    return;
    
}

void PairCopulaInvVfun_VecPar(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaInvVfun_VecPar(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaInvHfun_VecPar(family, theta, V, U, u, n);
            break;
        }
        case 180:
        {
            PairCopulaInvVfun_VecPar(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaInvHfun_VecPar(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
    }
    
    return;
    
}

double PairCopulaHfun(int family, const double theta, double U, double V)
{
    double u;
    PairCopulaHfun(family, &theta, &U, &V, &u, 1);
    return u;
}
double PairCopulaHfun(int family, const double theta0, const double theta1, double U, double V)
{
    double u;
    std::vector<double> theta(2);
    theta[0] = theta0;
    theta[1] = theta1;
    PairCopulaHfun(family, &theta[0], &U, &V, &u, 1);
    return u;
}
double PairCopulaHfun(int family, const double theta0, const double theta1, const double theta2, double U, double V)
{
    double u;
    std::vector<double> theta(3);
    theta[0] = theta0;
    theta[1] = theta1;
    theta[2] = theta2;
    PairCopulaHfun(family, &theta[0], &U, &V, &u, 1);
    return u;
}

double PairCopulaVfun(int family, const double theta, double U, double V)
{
    double u;
    PairCopulaVfun(family, &theta, &U, &V, &u, 1);
    return u;
}
double PairCopulaVfun(int family, const double theta0, const double theta1, double U, double V)
{
    double u;
    std::vector<double> theta(2);
    theta[0] = theta0;
    theta[1] = theta1;
    PairCopulaVfun(family, &theta[0], &U, &V, &u, 1);
    return u;
}
double PairCopulaVfun(int family, const double theta0, const double theta1, const double theta2, double U, double V)
{
    double u;
    std::vector<double> theta(3);
    theta[0] = theta0;
    theta[1] = theta1;
    theta[2] = theta2;
    PairCopulaVfun(family, &theta[0], &U, &V, &u, 1);
    return u;
}

// One Parameter
template <class T, class Policy>
        struct invhfun_bi
{
    invhfun_bi(T _family, T _theta, T _U, T _V)
    : family(_family), theta(_theta), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaHfun(family,theta,x,V) - U;
    }
        private:
            T family, theta, U, V;
};
// Two Parameters
template <class T, class Policy>
        struct invhfun_bi2
{
    invhfun_bi2(T _family, T _theta0, T _theta1, T _U, T _V)
    : family(_family), theta0(_theta0), theta1(_theta1), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaHfun(family,theta0,theta1,x,V) - U;
    }
        private:
            T family, theta0, theta1, theta2, U, V;
};
// Three Parameters
template <class T, class Policy>
        struct invhfun_bi3
{
    invhfun_bi3(T _family, T _theta0, T _theta1, T _theta2, T _U, T _V)
    : family(_family), theta0(_theta0), theta1(_theta1), theta2(_theta2), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaHfun(family,theta0,theta1,theta2,x,V) - U;
    }
        private:
            T family, theta0, theta1, theta2, U, V;
};

// One Parameter
template <class T, class Policy>
        struct invvfun_bi
{
    invvfun_bi(T _family, T _theta, T _U, T _V)
    : family(_family), theta(_theta), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaVfun(family,theta,U,x) - V;
    }
        private:
            T family, theta, U, V;
};
// Two Parameters
template <class T, class Policy>
        struct invvfun_bi2
{
    invvfun_bi2(T _family, T _theta0, T _theta1, T _U, T _V)
    : family(_family), theta0(_theta0), theta1(_theta1), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaVfun(family,theta0,theta1,U,x) - V;
    }
        private:
            T family, theta0, theta1, theta2, U, V;
};
// Three Parameters
template <class T, class Policy>
        struct invvfun_bi3
{
    invvfun_bi3(T _family, T _theta0, T _theta1, T _theta2, T _U, T _V)
    : family(_family), theta0(_theta0), theta1(_theta1), theta2(_theta2), U(_U), V(_V) {}
    
    T operator()(const T& x)
	{
        return PairCopulaVfun(family,theta0,theta1,theta2,U,x) - V;
    }
        private:
            T family, theta0, theta1, theta2, U, V;
};

double invhfun_bisect(int family, const double theta, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invhfun_bi<double, pol>(family,theta,U,V), min, max, tol).first;
}
double invhfun_bisect(int family, const double theta0, const double theta1, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invhfun_bi2<double, pol>(family,theta0,theta1,U,V), min, max, tol).first;
}
double invhfun_bisect(int family, const double theta0, const double theta1, const double theta2, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invhfun_bi3<double, pol>(family,theta0,theta1,theta2,U,V), min, max, tol).first;
}

double invvfun_bisect(int family, const double theta, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invvfun_bi<double, pol>(family,theta,U,V), min, max, tol).first;
}
double invvfun_bisect(int family, const double theta0, const double theta1, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invvfun_bi2<double, pol>(family,theta0,theta1,U,V), min, max, tol).first;
}
double invvfun_bisect(int family, const double theta0, const double theta1, const double theta2, double U, double V, unsigned int n)
{
    typedef boost::math::policies::policy<> pol;
    int digits = std::numeric_limits<double>::digits;
    boost::math::tools::eps_tolerance<double> tol(digits);
    
    double min = 1e-10;
    double max = 1-1e-10;
    return boost::math::tools::bisect(invvfun_bi3<double, pol>(family,theta0,theta1,theta2,U,V), min, max, tol).first;
}

void PairCopulaInvHfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 1: case 2: case 8: case 11: case 13: case 14: case 15:
        {
            // 1 Parameter: AMH, AsymFGM, FGM, Gumbel, Joe, PartialFrank, Plackett
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,*theta,UU,VV,n);
                
            }
            break;
        }
        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
        {
            // 2 Parameters: BB1, BB6, BB7, BB8, IteratedFGM, Tawn1, Tawn2
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,theta[0],theta[1],UU,VV,n);
                
            }
            break;
        }
        case 18:
        {
            // 3 Parameters: Tawn
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,theta[0],theta[1],theta[2],UU,VV,n);
                
            }
            break;
        }
        case 0:
        {
            // Indep
            for (i=0;i<n;i++)
            {
                *(u+i) = U[i];
            }
            break;
        }
        case 7:
        {
            //Clayton
            double h1,h2,UU,VV,uu;
            
            h1 = -1/ *theta;
            h2 = -*theta/(1+*theta);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = pow(pow(VV,-*theta)*(pow(UU,h2)-1)+1,h1);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 9:
        {
            // Frank
            double h1,h2,h3,h4,VV,UU,uu;
            h1 = exp(-*theta);
            h2 = expm1(-*theta);
            h3 = -1/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h4 = pow(h1,VV);
                
                uu = h3*log(1+h2/(h4*(1/UU-1)+1));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 10:
        {
            // Gaussian
            double x, y, h1,UU,VV,uu;
            
            h1 = sqrt(1-pow(*theta,2));
            
            boost::math::normal dist(0,1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                x = boost::math::quantile(dist, UU);
                y = boost::math::quantile(dist, VV);
                
                uu = boost::math::cdf(dist, x*h1+*theta*y);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 19:
        {
            // t
            double x, y, h1, nu1,UU,VV,uu;
            
            h1 = 1-pow(theta[0],2);
            nu1 = theta[1]+1;
            
            boost::math::students_t dist1(theta[1]);
            boost::math::students_t dist2(nu1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                x = boost::math::quantile(dist2, UU);
                y = boost::math::quantile(dist1, VV);
                
                uu = boost::math::cdf(dist1, x*sqrt((theta[1]+pow(y,2))*h1/nu1)+theta[0]*y);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
    }
    
    return;
}

void PairCopulaInvVfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 2:
        {
            // 1 Parameter: AsymFGM
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,*theta,UU,VV,n);
                
            }
            break;
        }
        case 16: case 17:
        {
            // 2 Parameters: Tawn1, Tawn2
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,theta[0],theta[1],UU,VV,n);
                
            }
            break;
        }
        case 18:
        {
            // 3 Parameters: Tawn
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,theta[0],theta[1],theta[2],UU,VV,n);
                
            }
            break;
        }
        default:
        {
            // symmetric copulas
            PairCopulaInvHfun(family, theta, V, U, u, n);
            break;
        }
    }
    return;
}


void PairCopulaInvHfun_VecPar(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 1: case 2: case 8: case 11: case 13: case 14: case 15:
        {
            // 1 Parameter: AMH, AsymFGM, FGM, Gumbel, Joe, PartialFrank, Plackett
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,theta[i],UU,VV,n);
                
            }
            break;
        }
        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
        {
            // 2 Parameters: BB1, BB6, BB7, BB8, IteratedFGM, Tawn1, Tawn2
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,theta[i],theta[i+n],UU,VV,n);
                
            }
            break;
        }
        case 18:
        {
            // 3 Parameters: Tawn
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invhfun_bisect(family,theta[i],theta[i+n],theta[i+2*n],UU,VV,n);
                
            }
            break;
        }
        case 0:
        {
            // Indep
            for (i=0;i<n;i++)
            {
                *(u+i) = U[i];
            }
            break;
        }
        case 7:
        {
            //Clayton
            double h1,h2,UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                h1 = -1/ theta[i];
                h2 = -theta[i]/(1+theta[i]);
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = pow(pow(VV,-theta[i])*(pow(UU,h2)-1)+1,h1);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 9:
        {
            // Frank
            double h1,h2,h3,h4,VV,UU,uu;
            
            for (i=0;i<n;i++)
            {
                h1 = exp(-theta[i]);
                h2 = expm1(-theta[i]);
                h3 = -1/ theta[i];
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h4 = pow(h1,VV);
                
                uu = h3*log(1+h2/(h4*(1/UU-1)+1));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 10:
        {
            // Gaussian
            double x, y, h1,UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                h1 = sqrt(1-pow(theta[i],2));
                
                boost::math::normal dist(0,1);
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                x = boost::math::quantile(dist, UU);
                y = boost::math::quantile(dist, VV);
                
                uu = boost::math::cdf(dist, x*h1+theta[i]*y);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 19:
        {
            // t
            double x, y, h1, nu1,UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                h1 = 1-pow(theta[i],2);
                nu1 = theta[i+n]+1;
                
                boost::math::students_t dist1(theta[i+n]);
                boost::math::students_t dist2(nu1);
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                x = boost::math::quantile(dist2, UU);
                y = boost::math::quantile(dist1, VV);
                
                uu = boost::math::cdf(dist1, x*sqrt((theta[i+n]+pow(y,2))*h1/nu1)+theta[i]*y);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
    }
    
    return;
}

void PairCopulaInvVfun_VecPar(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 2:
        {
            // 1 Parameter: AsymFGM
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,theta[i],UU,VV,n);
                
            }
            break;
        }
        case 16: case 17:
        {
            // 2 Parameters: Tawn1, Tawn2
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,theta[i],theta[i+n],UU,VV,n);
                
            }
            break;
        }
        case 18:
        {
            // 3 Parameters: Tawn
            
            double UU, VV;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u[i] = invvfun_bisect(family,theta[i],theta[i+n],theta[i+2*n],UU,VV,n);
                
            }
            break;
        }
        default:
        {
            // symmetric copulas
            PairCopulaInvHfun(family, theta, V, U, u, n);
            break;
        }
    }
    return;
}
