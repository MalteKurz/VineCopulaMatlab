#include "VineCPP_header.hpp"

double VineCopulaNegLL(VineCopula* Vine, double *U, int CutOffTree, unsigned int n)
{
    int i,j;
    double CLL=0;
    
    int d = Vine->Dimension;
    int J = min(CutOffTree,d-1);
    
    if (Vine->Type == 0)
    {
        std::vector<double> V((d-1)*n);
        if (Vine->Rotations==NULL)
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
                CLL += CLL1;
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=d-j-2;i>=0;i--)
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
                    CLL += CLL1;
                }
                for (i=d-j-2;i>=0;i--)
                {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], &V[(j+i)*n], n);
                }
            }
        }
        else
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[0],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                    CLL += CLL1;
                    PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0] , &V[i*n], n);
                }
                else
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
                    CLL += CLL1;
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=d-j-2;i>=0;i--)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&V[(j-1)*n],&V[(j+i)*n],&U1[0],&V1[0],rotation,n);
                        double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                        CLL += CLL1;
                    }
                    else
                    {
                        double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
                        CLL += CLL1;
                    }
                }
                for (i=d-j-2;i>=0;i--)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&V[(j-1)*n],&V[(j+i)*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(j+i)*n], n);
                    }
                    else
                    {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], &V[(j+i)*n], n);
                    }
                }
            }
        }
    }
    else
    {
        std::vector<double> V((d-2)*n);
        std::vector<double> H((d-2)*n);
        if (Vine->Rotations==NULL)
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
                CLL += CLL1;
                if (i<d-2) {
                    PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                }
                if (i>0) {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=0;i<d-j-1;i++)
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
                    CLL += CLL1;
                }
                for (i=0;i<d-j-1;i++)
                {
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &V[(i-1)*n], n);
                    }
                    if (i<d-j-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &H[i*n], n);
                    }
                }
            }
        }
        else
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[i*n],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                    CLL += CLL1;
                    if (i<d-2) {
                        PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                    }
                }
                else
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
                    CLL += CLL1;
                    if (i<d-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                    }
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=0;i<d-j-1;i++)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotation,n);
                        double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                        CLL += CLL1;
                    }
                    else
                    {
                        double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
                        CLL += CLL1;
                    }
                }
                for (i=0;i<d-j-1;i++)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotation,n);
                        if (i<d-j-2) {
                            PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
                        }
                        if (i>0) {
                            PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                        }
                    }
                    else
                    {
                        if (i>0) {
                            PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &V[(i-1)*n], n);
                        }
                        if (i<d-j-2) {
                            PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &H[i*n], n);
                        }
                    }
                }
            }
        }
    }
    
    return CLL;
}


double VineCopulaNegLL(const double *Thetas,VineCopula* Vine, double *U, int CutOffTree, unsigned int n)
{
    int i,j;
    double CLL=0;
    
    int d = Vine->Dimension;
    int J = min(CutOffTree,d-1);
    
    if (Vine->Type == 0)
    {
        std::vector<double> V((d-1)*n);
        if (Vine->Rotations==NULL)
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
                CLL += CLL1;
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=d-j-2;i>=0;i--)
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
                    CLL += CLL1;
                }
                for (i=d-j-2;i>=0;i--)
                {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], &V[(j+i)*n], n);
                }
            }
        }
        else
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[0],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                    CLL += CLL1;
                    PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0] , &V[i*n], n);
                }
                else
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
                    CLL += CLL1;
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=d-j-2;i>=0;i--)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&V[(j-1)*n],&V[(j+i)*n],&U1[0],&V1[0],rotation,n);
                        double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                        CLL += CLL1;
                    }
                    else
                    {
                        double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
                        CLL += CLL1;
                    }
                }
                for (i=d-j-2;i>=0;i--)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&V[(j-1)*n],&V[(j+i)*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(j+i)*n], n);
                    }
                    else
                    {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], &V[(j+i)*n], n);
                    }
                }
            }
        }
    }
    else
    {
        std::vector<double> V((d-2)*n);
        std::vector<double> H((d-2)*n);
        if (Vine->Rotations==NULL)
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
                CLL += CLL1;
                if (i<d-2) {
                    PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                }
                if (i>0) {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=0;i<d-j-1;i++)
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
                    CLL += CLL1;
                }
                for (i=0;i<d-j-1;i++)
                {
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &V[(i-1)*n], n);
                    }
                    if (i<d-j-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &H[i*n], n);
                    }
                }
            }
        }
        else
        {
            #pragma omp parallel for private(i) reduction(+:CLL)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[i*n],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                    CLL += CLL1;
                    if (i<d-2) {
                        PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                    }
                }
                else
                {
                    double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
                    CLL += CLL1;
                    if (i<d-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                    }
                }
            }
            
            for (j=1;j<J;j++)
            {
                #pragma omp parallel for private(i) reduction(+:CLL)
                for (i=0;i<d-j-1;i++)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotation,n);
                        double CLL1 = PairCopulaNegLL_Rotated_Obs(Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], n);
                        CLL += CLL1;
                    }
                    else
                    {
                        double CLL1 = PairCopulaNegLL(Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
                        CLL += CLL1;
                    }
                }
                for (i=0;i<d-j-1;i++)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotation,n);
                        if (i<d-j-2) {
                            PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
                        }
                        if (i>0) {
                            PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                        }
                    }
                    else
                    {
                        if (i>0) {
                            PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &V[(i-1)*n], n);
                        }
                        if (i<d-j-2) {
                            PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], &H[i*n], n);
                        }
                    }
                }
            }
        }
    }
    
    return CLL;
}