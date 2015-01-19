#include "VineCPP_header.hpp"

void VineCopulaFitSeq(VineCopula* Vine, double *CLL, double *U, unsigned int n)
{
    int i,j=0;
    
    int d = Vine->Dimension;
    
    if (Vine->Type == 0)
    {
        std::vector<double> V((d-1)*n);
        if (Vine->Rotations==NULL)
        {
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
            }
            for (i=0;i<d-1;i++)
            {
                PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
            }
            
            for (j=1;j<d-1;j++)
            {
                #pragma omp parallel for private(i)
                for (i=d-j-2;i>=0;i--)
                {
                    PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
                }
                for (i=d-j-2;i>=0;i--)
                {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], &V[(j+i)*n], n);
                }
            }
        }
        else
        {
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[0],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaFit_Rotated_Obs(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], rotation, &U1[0], &V1[0], n);
                }
                else
                {
                    PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], n);
                }
            }
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[0],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0] , &V[i*n], n);
                }
                else
                {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[0], &U[(i+1)*n], &V[i*n], n);
                }
            }
            
            for (j=1;j<d-1;j++)
            {
                #pragma omp parallel for private(i)
                for (i=d-j-2;i>=0;i--)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;   
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&V[(j-1)*n],&V[(j+i)*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaFit_Rotated_Obs(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], rotation, &U1[0], &V1[0], n);
                    }
                    else
                    {
                        PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &V[(j-1)*n], &V[(j+i)*n], n);
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
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
            }
            for (i=0;i<d-1;i++)
            {
                if (i<d-2) {
                    PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                }
                if (i>0) {
                    PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                }
            }
            
            for (j=1;j<d-1;j++)
            {
                #pragma omp parallel for private(i)
                for (i=0;i<d-j-1;i++)
                {
                    PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
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
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[i*n],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    PairCopulaFit_Rotated_Obs(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], rotation, &U1[0], &V1[0], n);
                }
                else
                {
                    PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], n);
                }
            }
            for (i=0;i<d-1;i++)
            {
                if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                {
                    int rotation;
                    std::vector<double> U1(n),V1(n);
                    rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                    
                    Rotate_Obs(&U[i*n],&U[(i+1)*n],&U1[0],&V1[0],rotation,n);
                    if (i<d-2) {
                        PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                    }
                }
                else
                {
                    if (i<d-2) {
                        PairCopulaHfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun((int) Vine->Families[d*j-j*(j+1)/2+i], Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                    }
                }
            }
            
            for (j=1;j<d-1;j++)
            {
                #pragma omp parallel for private(i)
                for (i=0;i<d-j-1;i++)
                {
                    if (Vine->Rotations[d*j-j*(j+1)/2+i]>0)
                    {
                        int rotation;
                        std::vector<double> U1(n),V1(n);
                        rotation = (int) Vine->Rotations[d*j-j*(j+1)/2+i];
                        
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotation,n);
                        PairCopulaFit_Rotated_Obs(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], rotation, &U1[0], &V1[0], n);
                    }
                    else
                    {
                        PairCopulaFit(Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], Vine->Families[d*j-j*(j+1)/2+i], &H[i*n], &V[i*n], n);
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
                        if (i>0) {
                            PairCopulaVfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &V[(i-1)*n], n);
                        }
                        if (i<d-j-2) {
                            PairCopulaHfun_Rotated_Obs((int) Vine->Families[d*j-j*(j+1)/2+i], rotation, Vine->Thetas+Vine->NumbParams[d*j-j*(j+1)/2+i], &U1[0], &V1[0], &H[i*n], n);
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
    
    *CLL = VineCopulaNegLL(Vine, U, d-1, n);
    return;
}

typedef struct {
    VineCopula *Vine;
    double *U;
    unsigned int n;
} my_data;

double ObjectiveFunctionVine(unsigned n, const double *x, double *grad, void *data)
{
    my_data *data_for_func = (my_data *) data;
    double CLL = VineCopulaNegLL(x,data_for_func->Vine,data_for_func->U,data_for_func->Vine->Dimension-1,data_for_func->n);
    return CLL;
}

void VineCopulaFit(VineCopula* Vine, double *CLL, double *x0, double *U, unsigned int n)
{
    VineCopulaFitSeq(Vine,&CLL[0],U,n);
    
    std::vector<double> bounds(120);
    LoadBounds(&bounds[0]);
    
    nlopt::opt opt;
    
    int d = Vine->Dimension;
    int J=Vine->NumbParams[(d-1)*d/2];
    int i;
    
    std::vector<double> lb(J);
    std::vector<double> ub(J);
    std::vector<double> x(J);
    opt = nlopt::opt(nlopt::LN_BOBYQA,J);
    
    for (i=0;i<(d-1)*d/2;i++)
    {
        switch(Vine->NumbParams[i+1]-Vine->NumbParams[i]){
            case 0:
            {
                break;
            }
            case 3:
            {
                lb[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6];
                lb[Vine->NumbParams[i]+1] = bounds[Vine->Families[i]*6+2];
                lb[Vine->NumbParams[i]+2] = bounds[Vine->Families[i]*6+4];
                ub[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6+1];
                ub[Vine->NumbParams[i]+1] = bounds[Vine->Families[i]*6+3];
                ub[Vine->NumbParams[i]+2] = bounds[Vine->Families[i]*6+5];
                x[Vine->NumbParams[i]] = *(Vine->Thetas+Vine->NumbParams[i]);
                x0[Vine->NumbParams[i]] = x[Vine->NumbParams[i]];
                x[Vine->NumbParams[i]+1] = *(Vine->Thetas+Vine->NumbParams[i]+1);
                x0[Vine->NumbParams[i]+1] = x[Vine->NumbParams[i]+1];
                x[Vine->NumbParams[i]+2] = *(Vine->Thetas+Vine->NumbParams[i]+2);
                x0[Vine->NumbParams[i]+2] = x[Vine->NumbParams[i]+2];
                break;
            }
            case 2:
            {
                lb[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6];
                lb[Vine->NumbParams[i]+1] = bounds[Vine->Families[i]*6+2];
                ub[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6+1];
                ub[Vine->NumbParams[i]+1] = bounds[Vine->Families[i]*6+3];
                x[Vine->NumbParams[i]] = *(Vine->Thetas+Vine->NumbParams[i]);
                x0[Vine->NumbParams[i]] = x[Vine->NumbParams[i]];
                x[Vine->NumbParams[i]+1] = *(Vine->Thetas+Vine->NumbParams[i]+1);
                x0[Vine->NumbParams[i]+1] = x[Vine->NumbParams[i]+1];
                break;
            }
            default:
            {
                lb[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6];
                ub[Vine->NumbParams[i]] = bounds[Vine->Families[i]*6+1];
                x[Vine->NumbParams[i]] = *(Vine->Thetas+Vine->NumbParams[i]);
                x0[Vine->NumbParams[i]] = x[Vine->NumbParams[i]];
            }
            
        }
    }
    
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
                
    my_data data = {Vine,U,n};
    
    opt.set_min_objective(ObjectiveFunctionVine,&data);
    
    opt.set_xtol_rel(1e-5);
    
    double minf;
    //nlopt::result result = opt.optimize(x, minf);
    opt.optimize(x, minf);
    
    //mexPrintf("result code = %i \n", opt.last_optimize_result());
    
    void nlopt_destroy(nlopt_opt opt);
    
    for (i=0;i<J;i++)
    {
        Vine->Thetas[i] = x[i];
    }
    
    CLL[1] = minf;
    
    return;
}
