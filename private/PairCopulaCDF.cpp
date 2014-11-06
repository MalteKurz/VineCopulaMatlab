#include "VineCPP_header.hpp"

void PairCopulaCDF(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaCDF(family, theta, U, V, p, n);
            break;
        }
        case 90:
        {
            PairCopulaCDF(family, theta, V, U, p, n);
            for (i=0;i<n;i++) p[i]=V[i]-p[i];
            break;
        }
        case 180:
        {
            PairCopulaCDF(family, theta, U, V, p, n);
            for (i=0;i<n;i++) p[i]=1-U[i]-V[i]+p[i];
            break;
        }
        case 270:
        {
            PairCopulaCDF(family, theta, V, U, p, n);
            for (i=0;i<n;i++) p[i]=U[i]-p[i];
            break;
        }
    }
    
    return;
    
}

void PairCopulaCDF(int family, const double *theta, double *U, double *V, double *p, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 0:
        {
            // Indep
            for (i=0;i<n;i++) p[i]=U[i]*V[i];
            break;
        }
        case 1:
        {
            // AMH
            double h1,UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h1 = UU*VV;
                
                p[i] = h1/(1-*theta*(1-UU-VV+h1));
            }
            break;
        }
        case 2:
        {
            // AsymFGM
            double h1,UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h1 = UU*VV;
                
                p[i] = h1+*theta*h1*(1-h1)*VV*(1-UU);
            }
            break;
        }
        case 3:
        {
            // BB1
            double h1,h2,h3;
            
            h1 = -theta[0];
            h2 = 1/h1;
            h3 = 1/ theta[1];
            
            for (i=0;i<n;i++)
            {
                double hu,hv,huv,hhuv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(pow(UU,h1)-1,theta[1]);
                hv = pow(pow(VV,h1)-1,theta[1]);
                huv = hu + hv;
                hhuv = 1+pow(huv,h3);
                
                p[i] = pow(hhuv,h2);
            }
            break;
        }
        case 4:
        {
            // BB6
            double h1,h2,h3,h4,h5,h6;
            
            h1 = 1/theta[0];
            h2 = 1/theta[1];
            
            for (i=0;i<n;i++)
            {
                double u1,v1,ut,vt,hu,hv,huv,hhuv,ehhuv,hhhuv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u1 = 1-UU;
                v1 = 1-VV;
                ut = 1-pow(u1,theta[0]);
                vt = 1-pow(v1,theta[0]);
                hu = -log(ut);
                hv = -log(vt);
                huv = pow(hu,theta[1])+pow(hv,theta[1]);
                hhuv = pow(huv,h2);
                ehhuv = exp(-hhuv);
                hhhuv = 1-ehhuv;
                
                p[i] = 1-pow(hhhuv,h1);
            }
            break;
        }
        case 5:
        {
            // BB7
            double h1,h2,h3;
            
            h1 = -theta[1];
            h2 = 1/theta[0];
            h3 = 1/theta[1];
            
            for (i=0;i<n;i++)
            {
                double hu,hv,huv,htilde,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = 1-UU;
                hv = 1-VV;
                huv = pow(1-pow(hu,theta[0]),h1)+pow(1-pow(hv,theta[0]),h1)-1;
                htilde = pow(huv,h3);
                
                p[i] = 1-pow(1-1/htilde,h2);
            }
            break;
        }
        case 6:
        {
            // BB8
                double h1,h2,h3,h4,h5;
                
                h1 = 1-theta[1];
                h2 = pow(h1,theta[0]);
                h3 = 1/(1-h2);
                h4 = 1/theta[0];
                h5 = 1/theta[1];
                
                for (i=0;i<n;i++)
                {
                    double hu,hv,UU,VV;
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    hu = 1-pow(1-theta[1]*UU,theta[0]);
                    hv = 1-pow(1-theta[1]*VV,theta[0]);
                    
                    p[i] = h5*(1-pow(1-h3*hu*hv,h4));
                }
            break;
        }
        case 7:
        {
            //Clayton
            double h1,h2,hu,hv,UU,VV;
            
            h1 = -*theta;
            h2 = 1/h1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                hu = pow(UU,h1);
                hv = pow(VV,h1);
                
                p[i] = pow(hu+hv-1,h2);
            }
            break;
        }
        case 8:
        {
            // FGM
            double huv,UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                huv = UU*VV;
                
                p[i] = huv + *theta*huv*(1-huv);
            }
            break;
        }
        case 9:
        {
            // Frank
            double h1,h2,h3,h4,UU,VV;
            
            h1 = -*theta;
            h2 = expm1(h1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h3 = expm1(h1*UU)*expm1(h1*VV);
                
                p[i] = -log(1+h3/h2)/ *theta;
            }
            break;
        }
        case 10:
        {
            // Gaussian
            double UU,VV,CORREL,ABSEPS=0.001,RELEPS=0,ERROR=0,rho;
            int N=2,NU=0,MAXPTS=25000,INFORM;
            std::vector<double> ZeroArray(2), UPPER(2);
            std::vector<int> ZeroIntArray(2);
            
            rho=*theta;
            
            boost::math::normal dist(0,1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                UPPER[0] = boost::math::quantile(dist, UU);
                UPPER[1] = boost::math::quantile(dist, VV);
                mvtdst_(&N,&NU,&ZeroArray[0],&UPPER[0],&ZeroIntArray[0],&rho,&ZeroArray[0],&MAXPTS,&ABSEPS,&RELEPS,&ERROR,&p[i],&INFORM);
            }

            break;
        }
        case 11:
        {
            // Gumbel
            double h1,h2,h3,h4,h5,h7,UU,VV;
            
            h1 = 1/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h2 = -log(UU);
                h3 = -log(VV);
                h4 = pow(h2,*theta) + pow(h3,*theta);
                h5 = pow(h4,h1);
                
                p[i] = exp(-h5);
            }
            break;
        }
        case 12:
        {
            // IteratedFGM
            double pUV,p1UV,UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                pUV = UU*VV;
                p1UV = (1-UU)*(1-VV);
                
                p[i] = pUV*(1 + *theta*p1UV+theta[1]*pUV*p1UV);
            }
            break;
        }
        case 13:
        {
            // Joe
            double h1,hu,hv,huv,u1,v1,UU,VV;
            
            h1 = 1/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u1 = 1-UU;
                hu = pow(u1,*theta);
                v1 = 1-VV;
                hv = pow(v1,*theta);
                huv = hu+hv-hu*hv;
                
                p[i] = 1-pow(huv,h1);
            }
            break;
        }
        case 14:
        {
            // PartialFrank
            double h1,pUV,sUV,UU,VV;
            
            h1 = expm1(-*theta);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                pUV = UU*VV;
                sUV = UU+VV;
                
                p[i] = pUV/(*theta*(sUV-pUV))*(log(1+h1*(1-sUV+pUV))+*theta);
            }
            break;
        }
        case 15:
        {
            // Plackett
            if (*theta == 1)
            {
                for (i=0;i<n;i++)
                {
                    p[i] = U[i]*V[i];
                }
            }
            else
            {
                double h1,h2,h3,sUV,pUV,UU,VV;
                
                h1 = *theta-1;
                h2 = 4**theta*h1;
                h3 = 2*h1;
                
                for (i=0;i<n;i++)
                {
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    sUV = UU+VV;
                    pUV = UU*VV;
                    
                    p[i] = ((1+h1*sUV)-sqrt(pow(1+h1*sUV,2)-h2*pUV))/h3;
                }
            }
            break;
        }
        case 16:
        {
            // Tawn1
            std::vector<double> U1(n);
            for (i=0;i<n;i++) U1[i]=pow(U[i],theta[1]);
            
            PairCopulaCDF(11,theta,&U1[0],V,p,n);
            double h1;
            
            h1 = 1-theta[1];
            
            for (i=0;i<n;i++)
            {
                double hu,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,h1);
                
                p[i] *= hu;
            }
            break;
        }
        case 17:
        {
            // Tawn2
            std::vector<double> V1(n);
            for (i=0;i<n;i++) V1[i]=pow(V[i],theta[1]);
            
            PairCopulaCDF(11,theta,U,&V1[0],p,n);
            double h1;
            
            h1 = 1-theta[1];
            
            for (i=0;i<n;i++)
            {
                double hv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hv = pow(VV,h1);
                
                p[i] *= hv;
            }
            break;
        }
        case 18:
        {
            // Tawn
            std::vector<double> U1(n);
            std::vector<double> V1(n);
            for (i=0;i<n;i++) 
            {
                U1[i]=pow(U[i],theta[1]);
                V1[i]=pow(V[i],theta[2]);
            }
            
            PairCopulaCDF(11,theta,&U1[0],&V1[0],p,n);
            double h1,h2;
            
            h1 = 1-theta[1];
            h2 = 1-theta[2];
            
            for (i=0;i<n;i++)
            {
                double hu,hv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,h1);
                hv = pow(VV,h2);
                
                p[i] *= hu*hv;
            }
            break;
        }
        case 19:
        {
            // t
            double UU,VV,CORREL,ABSEPS=0.001,RELEPS=0,ERROR=0,rho;
            int N=2,nu,MAXPTS=25000,INFORM;
            std::vector<double> ZeroArray(2), UPPER(2);
            std::vector<int> ZeroIntArray(2);
            
            rho = theta[0];
            nu = (int) theta[1];
            
            if (nu != theta[1])
            {
                mexErrMsgTxt( "The CDF of the t copula is only implemented for integer valued degrees of freedom. Call the build-in MATLAB function copulacdf instead.");
            }
            
            boost::math::students_t dist(theta[1]);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                UPPER[0] = boost::math::quantile(dist, UU);
                UPPER[1] = boost::math::quantile(dist, VV);
                mvtdst_(&N,&nu,&ZeroArray[0],&UPPER[0],&ZeroIntArray[0],&rho,&ZeroArray[0],&MAXPTS,&ABSEPS,&RELEPS,&ERROR,&p[i],&INFORM);
            }
            
            break;
        }
        
    }
    
    return;
}
