#include "VineCPP_header.hpp"

void VineCopulaRand(VineCopula* Vine, double *U, unsigned int n)
{
    int k,l,i,j;
    
    boost::mt19937 gen;
    // Load the state
    std::ifstream fi(PathSeed);
    fi>>gen;
    fi.close();
    
    boost::uniform_01 <> URAND;
    
    boost::variate_generator <boost::mt19937&, boost::uniform_01 <> > RAND(gen,URAND);
    
    int d = Vine->Dimension;
    
    std::vector<double> W(d*n);
    
    for (j=0;j<d;j++)
    {
        for (i=0;i<(int)n;i++)
        {
            W[j*n+i] = RAND();
            U[j*n+i] = W[j*n+i];
        }
    }
    
    // Save the state
    std::ofstream fo(PathSeed,
            std::ios_base::out);
    fo<<gen;
    fo.close();
    
    
    if (Vine->Type == 0)
    {
        std::vector<int> ind(d*(d-1)/2);
        
        ind[0] = 0;
        ind[1] = d-1;
        ind[2] = 1;
        
        j=3;
        for (i=3;i<d;i++)
        {
            ind[j] = ind[j-i+1]+d-i+1;
            for (k=j+1;k<j+i;k++)
            {
                ind[k] = ind[k-i]+1;
            }
            j = j+i;
        }
        
        if (Vine->Rotations==NULL)
        {
            PairCopulaInvVfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &W[0], &U[n] , &U[n], n);
            
            k=1;
            for (i=2;i<d;i++)
            {
                for (j=i-1;j>=0;j--)
                {
                    PairCopulaInvVfun((int) Vine->Families[ind[k]], Vine->Thetas+Vine->NumbParams[ind[k]], &W[j*n], &U[i*n] , &U[i*n], n);
                    k++;
                }
            }
        }
        else
        {
            int rotation;
            std::vector<double> U1(n),V1(n);
            if (Vine->Rotations[0]>0)
            {
                rotation = (int) Vine->Rotations[0];
                
                Rotate_Obs(&W[0],&U[n],&U1[0],&V1[0],rotation,n);
                PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[0], rotation, Vine->Thetas+Vine->NumbParams[0], &U1[0], &V1[0] , &U[n], n);
            }
            else
            {
                PairCopulaInvVfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &W[0], &U[n] , &U[n], n);
            }
            
            k=1;
            for (i=2;i<d;i++)
            {
                for (j=i-1;j>=0;j--)
                {
                    if (Vine->Rotations[ind[k]]>0)
                    {
                        rotation = (int) Vine->Rotations[ind[k]];
                        
                        Rotate_Obs(&W[j*n],&U[i*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[ind[k]], rotation, Vine->Thetas+Vine->NumbParams[ind[k]], &U1[0], &V1[0] , &U[i*n], n);
                    }
                    else
                    {
                        PairCopulaInvVfun((int) Vine->Families[ind[k]], Vine->Thetas+Vine->NumbParams[ind[k]], &W[j*n], &U[i*n] , &U[i*n], n);
                    }
                    k++;
                }
            }
        }
    }
    else
    {
        std::vector<int> indInvVfuns(d*(d-1)/2);
        std::vector<int> indHfuns(d*(d-1)/2);
        
        std::vector<double> A((d-2)*n);
        std::vector<double> B((d-2)*n);
        
        indInvVfuns[0] = 0;
        indInvVfuns[1] = d-1;
        indInvVfuns[2] = 1;
        indHfuns[0] = 0;
        indHfuns[1] = 1;
        indHfuns[2] = d-1;
        
        
        j=3;
        for (i=3;i<d;i++)
        {
            indInvVfuns[j] = indInvVfuns[j-i+1]+d-i+1;
            indHfuns[j+i-1] = indHfuns[j-1]+d-i+1;
            for (k=j+1;k<j+i;k++)
            {
                indInvVfuns[k] = indInvVfuns[k-i]+1;
                indHfuns[k-1] = indHfuns[k-i]+1;
            }
            j = j+i;
        }
        
        
        
        if (Vine->Rotations==NULL)
        {
            PairCopulaInvVfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &W[0], &W[n] , &U[n], n);
            PairCopulaHfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &U[0], &U[n] , &B[0], n);
            
            k=1,l=1;
            for (i=2;i<d-1;i++)
            {
                PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[0], &W[i*n] , &A[0], n);
                k++;
                for (j=2;j<i;j++)
                {
                    PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[(j-1)*n], &A[(j-2)*n] , &A[(j-1)*n], n);
                    k++;
                }
                PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U[(i-1)*n], &A[(i-2)*n] , &U[i*n], n);
                k++;
                
                
                PairCopulaHfun((int) Vine->Families[indHfuns[l]], Vine->Thetas+Vine->NumbParams[indHfuns[l]], &U[(i-1)*n], &U[i*n] , &B[(i-1)*n], n);
                l++;
                
                for (j=i-2;j>=0;j--)
                {
                    PairCopulaHfun((int) Vine->Families[indHfuns[l]], Vine->Thetas+Vine->NumbParams[indHfuns[l]], &B[j*n], &A[j*n] , &B[j*n], n);
                    l++;
                }
            }
            i=d-1;
            PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[0], &W[i*n] , &A[0], n);
            k++;
            for (j=2;j<i;j++)
            {
                PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[(j-1)*n], &A[(j-2)*n] , &A[(j-1)*n], n);
                k++;
            }
            PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U[(i-1)*n], &A[(i-2)*n] , &U[i*n], n);
            k++;
        }
        else
        {
            int rotation;
            std::vector<double> U1(n),V1(n);
            
            if (Vine->Rotations[0]>0)
            {
                rotation = (int) Vine->Rotations[0];
                
                Rotate_Obs(&W[0],&W[n],&U1[0],&V1[0],rotation,n);
                PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[0], rotation, Vine->Thetas+Vine->NumbParams[0], &U1[0], &V1[0] , &U[n], n);
                
                Rotate_Obs(&U[0],&U[n],&U1[0],&V1[0],rotation,n);
                PairCopulaHfun_Rotated_Obs((int) Vine->Families[0], rotation, Vine->Thetas+Vine->NumbParams[0], &U1[0], &V1[0] , &B[0], n);
            }
            else
            {
                PairCopulaInvVfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &W[0], &W[n] , &U[n], n);
                PairCopulaHfun((int) Vine->Families[0], Vine->Thetas+Vine->NumbParams[0], &U[0], &U[n] , &B[0], n);
            }
            
            k=1,l=1;
            for (i=2;i<d-1;i++)
            {
                if (Vine->Rotations[indInvVfuns[k]]>0)
                {
                    rotation = (int) Vine->Rotations[indInvVfuns[k]];
                    
                    Rotate_Obs(&B[0],&W[i*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &A[0], n);
                }
                else
                {
                    PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[0], &W[i*n] , &A[0], n);
                }
                k++;
                for (j=2;j<i;j++)
                {
                    if (Vine->Rotations[indInvVfuns[k]]>0)
                    {
                        rotation = (int) Vine->Rotations[indInvVfuns[k]];
                        
                        Rotate_Obs(&B[(j-1)*n],&A[(j-2)*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &A[(j-1)*n], n);
                    }
                    else
                    {
                        PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[(j-1)*n], &A[(j-2)*n] , &A[(j-1)*n], n);
                    }
                    k++;
                }
                if (Vine->Rotations[indInvVfuns[k]]>0)
                {
                    rotation = (int) Vine->Rotations[indInvVfuns[k]];
                    
                    Rotate_Obs(&U[(i-1)*n],&A[(i-2)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &U[i*n], n);
                }
                else
                {
                    PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U[(i-1)*n], &A[(i-2)*n] , &U[i*n], n);
                }
                k++;
                
                
                if (Vine->Rotations[indHfuns[l]]>0)
                {
                    rotation = (int) Vine->Rotations[indHfuns[l]];
                    
                    Rotate_Obs(&U[(i-1)*n],&U[i*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaHfun_Rotated_Obs((int) Vine->Families[indHfuns[l]], rotation, Vine->Thetas+Vine->NumbParams[indHfuns[l]], &U1[0], &V1[0] , &B[(i-1)*n], n);
                }
                else
                {
                    PairCopulaHfun((int) Vine->Families[indHfuns[l]], Vine->Thetas+Vine->NumbParams[indHfuns[l]], &U[(i-1)*n], &U[i*n] , &B[(i-1)*n], n);
                }
                l++;
                
                for (j=i-2;j>=0;j--)
                {
                    if (Vine->Rotations[indHfuns[l]]>0)
                    {
                        rotation = (int) Vine->Rotations[indHfuns[l]];
                        
                        Rotate_Obs(&B[j*n],&A[j*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaHfun_Rotated_Obs((int) Vine->Families[indHfuns[l]], rotation, Vine->Thetas+Vine->NumbParams[indHfuns[l]], &U1[0], &V1[0] , &B[j*n], n);
                    }
                    else
                    {
                        PairCopulaHfun((int) Vine->Families[indHfuns[l]], Vine->Thetas+Vine->NumbParams[indHfuns[l]], &B[j*n], &A[j*n] , &B[j*n], n);
                    }
                    l++;
                }
            }
            i=d-1;
            if (Vine->Rotations[indInvVfuns[k]]>0)
            {
                rotation = (int) Vine->Rotations[indInvVfuns[k]];
                
                Rotate_Obs(&B[0],&W[i*n],&U1[0],&V1[0],rotation,n);
                PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &A[0], n);
            }
            else
            {
                PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[0], &W[i*n] , &A[0], n);
            }
            k++;
            for (j=2;j<i;j++)
            {
                if (Vine->Rotations[indInvVfuns[k]]>0)
                {
                    rotation = (int) Vine->Rotations[indInvVfuns[k]];
                    
                    Rotate_Obs(&B[(j-1)*n],&A[(j-2)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &A[(j-1)*n], n);
                }
                else
                {
                    PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &B[(j-1)*n], &A[(j-2)*n] , &A[(j-1)*n], n);
                }
                k++;
            }
            if (Vine->Rotations[indInvVfuns[k]]>0)
            {
                rotation = (int) Vine->Rotations[indInvVfuns[k]];
                
                Rotate_Obs(&U[(i-1)*n],&A[(i-2)*n],&U1[0],&V1[0],rotation,n);
                PairCopulaInvVfun_Rotated_Obs((int) Vine->Families[indInvVfuns[k]], rotation, Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U1[0], &V1[0] , &U[i*n], n);
            }
            else
            {
                PairCopulaInvVfun((int) Vine->Families[indInvVfuns[k]], Vine->Thetas+Vine->NumbParams[indInvVfuns[k]], &U[(i-1)*n], &A[(i-2)*n] , &U[i*n], n);
            }
            k++;
        }
    }
    
    return;
}

