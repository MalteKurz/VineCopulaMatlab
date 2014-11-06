#include "VineCPP_header.hpp"

void PairCopulaSelect(int *family, double *theta, int *rotation, double *U, double *V, unsigned int n, double *familyset, int m)
{
    // Testing for the independence copula
    double *tau, AIC0=0,AIC1=0,Tau=0;
    int j;
    
    tau = &Tau;
    
    PairCopulaIndepTest(U,V,n,tau);
    
    if (PairCopulaIndepTest(U,V,n,tau))
    {
        *family = 0;
        *rotation = 0;
    }
    else
    {
        if (*tau<0)
        {
            std::vector<double> U90(n),V90(n),U270(n),V270(n);
            
            Rotate_Obs(U,V,&U90[0],&V90[0],90,n);
            Rotate_Obs(U,V,&U270[0],&V270[0],270,n);
            
            if (*familyset == -1 && m == 1)
            {
                #pragma omp parallel for private(j,AIC1)
                for (j=1;j<20;j++)
                {
                    switch(j){
                        case 8: case 9: case 10: case 15:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            break;
                        }
                        case 19:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        case 18:
                        {
                            std::vector<double> Theta(3);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            break;
                        }    
                        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }    
                        default:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                #pragma omp parallel for private(j,AIC1)
                for (j=0;j<m;j++)
                {
                    int J=familyset[j];
                    switch(J){
                        case 8: case 9: case 10: case 15:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            break;
                        }
                        case 19:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        case 18:
                        {
                            std::vector<double> Theta(3);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            break;
                        } 
                        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        default:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,90,&U90[0],&V90[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 90;
                                    theta[0] = Theta[0];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,270,&U270[0],&V270[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 270;
                                    theta[0] = Theta[0];
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            std::vector<double> U180(n),V180(n);
            
            Rotate_Obs(U,V,&U180[0],&V180[0],180,n);
            
            if (*familyset == -1 && m == 1)
            {
                #pragma omp parallel for private(j,AIC1)
                for (j=1;j<20;j++)
                {
                    switch(j){
                        case 8: case 9: case 10: case 15:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            break;
                        }
                        case 19:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        case 18:
                        {
                            std::vector<double> Theta(3);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            break;
                        }    
                        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        default:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],j,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],j,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = j;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                #pragma omp parallel for private(j,AIC1)
                for (j=0;j<m;j++)
                {
                    int J=familyset[j];
                    switch(J){
                        case 8: case 9: case 10: case 15:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            break;
                        }
                        case 19:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        case 18:
                        {
                            std::vector<double> Theta(3);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                    theta[2] = Theta[2];
                                }
                            }
                            break;
                        }
                        case 3: case 4: case 5: case 6: case 12: case 16: case 17:
                        {
                            std::vector<double> Theta(2);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                    theta[1] = Theta[1];
                                }
                            }
                            break;
                        }
                        default:
                        {
                            std::vector<double> Theta(1);
                            PairCopulaAIC(&AIC1,&Theta[0],J,U,V,n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 0;
                                    theta[0] = Theta[0];
                                }
                            }
                            
                            PairCopulaAIC_Rotated_Obs(&AIC1,&Theta[0],J,180,&U180[0],&V180[0],n);
                            #pragma omp critical(FamilyUpdate)
                            {
                                if (AIC1<AIC0)
                                {
                                    AIC0=AIC1;
                                    *family = J;
                                    *rotation = 180;
                                    theta[0] = Theta[0];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}
