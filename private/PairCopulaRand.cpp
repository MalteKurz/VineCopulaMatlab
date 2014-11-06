#include "VineCPP_header.hpp"
void PairCopulaRand(int family, int rotation, const double *theta, double *U, double *V, unsigned int n)
{     
    unsigned int i;
    
    boost::mt19937 gen;
    // Load the state
    std::ifstream fi(PathSeed);
    fi>>gen;
    fi.close();
    
    boost::uniform_01 <> URAND;
    
    boost::variate_generator <boost::mt19937&, boost::uniform_01 <> > RAND(gen,URAND);
    
    for (i=0;i<n;i++)
    {
        U[i] = RAND();
        V[i] = RAND();
    }
    
    // Save the state
    std::ofstream fo(PathSeed,
            std::ios_base::out);
    fo<<gen;
    fo.close();
    
    switch(family){
        case 0:
        {
            // Indep
            break;
        }
        default:
        {
            // AMH, AsymFGM, BB6, BB7, Gaussian, Gumbel, IteratedFGM, Joe, Plackett, Tawn1, Tawn2, Tawn, t
            if(rotation>0)
            {
                std::vector<double> U1(n),V1(n);
                
                Rotate_Obs(U,V,&U1[0],&V1[0],rotation,n);
                PairCopulaInvHfun_Rotated_Obs(family, rotation, theta, &U1[0], &V1[0], U, n);
            }
            else
            {
                PairCopulaInvHfun(family, rotation, theta, U, V, U, n);
            }
            break;
        }
        
    }
    
    return;
}

void PairCopulaRand(int family, const double *theta, double *U, double *V, unsigned int n)
{
    unsigned int i;
    
    boost::mt19937 gen;
    // Load the state
    std::ifstream fi(PathSeed);
    fi>>gen;
    fi.close();
    
    boost::uniform_01 <> URAND;
    
    boost::variate_generator <boost::mt19937&, boost::uniform_01 <> > RAND(gen,URAND);
    
    for (i=0;i<n;i++)
    {
        U[i] = RAND();
        V[i] = RAND();
    }
    
    // Save the state
    std::ofstream fo(PathSeed,
            std::ios_base::out);
    fo<<gen;
    fo.close();
    
    switch(family){
        case 0:
        {
            // Indep
            break;
        }
        default:
        {
            // AMH, AsymFGM, BB6, BB7, Gaussian, Gumbel, IteratedFGM, Joe, Plackett, Tawn1, Tawn2, Tawn, t
            PairCopulaInvHfun(family, theta, U, V, U, n);
            break;
        }
        
    }
    
    return;
}

