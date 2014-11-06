#include "VineCPP_header.hpp"

void VineCopulaGetPseudoObs(VineCopula* Vine, double *U, double *V, unsigned int n)
{
    int l,i,j=0;
    
    int d = Vine->Dimension;
    
    if (Vine->Rotations==NULL)
    {
        for (i=0;i<d-1;i++)
        {
            PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[(d*j-j*(j+1)/2+i)*n], n);
        }
        
        for (j=1;j<d-2;j++)
        {
            for (i=d-j-2;i>=0;i--)
            {
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(d*(j-1)-j*(j-1)/2)*n], &V[(d*(j-1)-j*(j-1)/2+i+1)*n], &V[(d*j-j*(j+1)/2+i)*n], n);
            }
        }
    }
    else
    {
        for (i=0;i<d-1;i++)
        {
            if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
            {
                int rotation;
                std::vector<double> U1(n),V1(n);
                rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                
                Rotate_Obs(&U[0],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(d*j-j*(j+1)/2+i)*n], n);
            }
            else
            {
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[(d*j-j*(j+1)/2+i)*n], n);
            }
        }
        
        for (j=1;j<d-2;j++)
        {
            for (i=d-j-2;i>=0;i--)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&V[(d*(j-1)-j*(j-1)/2)*n],&V[(d*(j-1)-j*(j-1)/2+i+1)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(d*j-j*(j+1)/2+i)*n], n);
                }
                else
                {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(d*(j-1)-j*(j-1)/2)*n], &V[(d*(j-1)-j*(j-1)/2+i+1)*n], &V[(d*j-j*(j+1)/2+i)*n], n);
                }
            }
        }
    }
    
    return;
}

void VineCopulaGetPseudoObs(VineCopula* Vine, double *U, double *H, double *V, unsigned int n)
{
    int l,i,j=0;
    
    int d = Vine->Dimension;
    
    if (Vine->Rotations==NULL)
    {
        for (i=0;i<d-1;i++)
        {
            if (i<d-2) {
                PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
            }
            if (i>0) {
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
            }
        }
        
        for (j=1;j<d-2;j++)
        {
            for (i=0;i<d-j-1;i++)
            {
                if (i>0) {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*j-j*(j+1)/2+i-j-1)*n], n);
                }
                if (i<d-j-2) {
                    PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &H[(d*j-j*(j+1)/2+i-j)*n], n);
                }
            }
        }
    }
    else
    {
        for (i=0;i<d-1;i++)
        {
            if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
            {
                int rotation;
                std::vector<double> U1(n),V1(n);
                rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                
                Rotate_Obs(&U[i*n], &U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                if (i<d-2) {
                    PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[(d*j-j*(j+1)/2+i-j)*n], n);
                }
                if (i>0) {
                    PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(d*j-j*(j+1)/2+i-j-1)*n], n);
                }
            }
            else
            {
                if (i<d-2) {
                    PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[(d*j-j*(j+1)/2+i-j)*n], n);
                }
                if (i>0) {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(d*j-j*(j+1)/2+i-j-1)*n], n);
                }
            }
        }
        
        for (j=1;j<d-2;j++)
        {
            for (i=0;i<d-j-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&H[(d*(j-1)-j*(j-1)/2+i-j+1)*n],&V[(d*(j-1)-j*(j-1)/2+i-j+1)*n],&U1[0],&V1[0],rotation,n);
                    if (i>0) {
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(d*j-j*(j+1)/2+i-j-1)*n], n);
                    }
                    if (i<d-j-2) {
                        PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[(d*j-j*(j+1)/2+i-j)*n], n);
                    }
                }
                else
                {
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*j-j*(j+1)/2+i-j-1)*n], n);
                    }
                    if (i<d-j-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &V[(d*(j-1)-j*(j-1)/2+i-j+1)*n], &H[(d*j-j*(j+1)/2+i-j)*n], n);
                    }
                }
            }
        }
    }
    
    return;
}
