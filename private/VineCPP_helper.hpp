#include "VineCPP_header.hpp"

#define LowerBoundUnitInterval  1e-10
#define UpperBoundUnitInterval  1-1e-10

inline double max(double x, double y)
{
    return (x > y)? x : y;
}

inline unsigned int max(unsigned int x, unsigned int y)
{
    return (x > y)? x : y;
}

inline double min(double x, double y)
{
    return (x > y)? y : x;
}

inline int IndTria(int i, int j, int n)
{
    return (j > i)? (i*n-(i-1)*i/2+j-i) : (j*n-(j-1)*j/2+i-j);
}

inline double sqrt( int val )
{
    return sqrt( (double) val );
}

inline double CheckBounds(double u)
{
    return max(min(u,UpperBoundUnitInterval),LowerBoundUnitInterval);
}

void Rotate_Obs(double *U,double *V, int rotation, unsigned int n);
void Rotate_Obs(double *U,double *V, double *U1, double *V1, int rotation, unsigned int n);

//void LogRand(double *X, double theta, unsigned int n);
//void StableRand(double *X, double alpha, double beta, double gamma, double delta, unsigned int n);

extern "C"{
    void mvtdst_(int *N,int *NU,double *LOWER,double *UPPER,int *INFIN,double *CORREL,double *DELTA,int *MAXPTS,double *ABSEPS,double *RELEPS,double *ERROR,double *VALUE,int *INFORM );
}

// Vine copula struct and its constructor
struct VineCopula{
    int Type; // The type (0=C-Vine, 1=D-Vine)
    int Dimension; // The dimension of the vine
    double *Families; // Pair-copula families
    double *Rotations; // Rotation degrees of the pair-copula families
    double *Thetas; // Pointer to the parameter vector
    std::vector<int> NumbParams; // Cumulative sum of the number of parameters of each pair-copula
    VineCopula(int type, int dimension, double *families, double *thetas)
    {
        Type = type;
        Dimension = dimension;
        Families = families;
        Thetas = thetas;
        Rotations = NULL;
        
        NumbParams.resize((dimension-1)*dimension/2+1);
            
        int i;
        for (i=0;i<(dimension-1)*dimension/2;i++)
        {
            switch((int) families[i]){
                case 0:
                {
                    NumbParams[i+1] = NumbParams[i];
                    break;
                }
                case 18:
                {
                    NumbParams[i+1] = NumbParams[i] +3;
                    break;
                }
                case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
                {
                    NumbParams[i+1] = NumbParams[i] +2;
                    break;
                }
                default:
                {
                    NumbParams[i+1] = NumbParams[i] +1;
                }
                
            }
        }
    }
    VineCopula(int type, int dimension, double *families, double *rotations, double *thetas)
    {
        Type = type;
        Dimension = dimension;
        Families = families;
        Thetas = thetas;
        Rotations = rotations;
        
        NumbParams.resize((dimension-1)*dimension/2+1);
            
        int i;
        for (i=0;i<(dimension-1)*dimension/2;i++)
        {
            switch((int) families[i]){
                case 0:
                {
                    NumbParams[i+1] = NumbParams[i];
                    break;
                }
                case 18:
                {
                    NumbParams[i+1] = NumbParams[i] +3;
                    break;
                }
                case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
                {
                    NumbParams[i+1] = NumbParams[i] +2;
                    break;
                }
                default:
                {
                    NumbParams[i+1] = NumbParams[i] +1;
                }
                
            }
        }
    }
};
