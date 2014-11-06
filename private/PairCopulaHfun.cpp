#include "VineCPP_header.hpp"

void PairCopulaHfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{  
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaHfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaVfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 180:
        {
            PairCopulaHfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaVfun(family, theta, V, U, u, n);
            break;
        }
    }
    
    return;
    
}

void PairCopulaVfun_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaVfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaHfun(family, theta, V, U, u, n);
            break;
        }
        case 180:
        {
            PairCopulaVfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaHfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
    }
    
    return;
    
}

void PairCopulaHfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaHfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaVfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 180:
        {
            PairCopulaHfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaVfun(family, theta, V, U, u, n);
            break;
        }
    }
    
    return;
    
}

void PairCopulaVfun(int family, int rotation, const double *theta, double *U, double *V, double *u, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    unsigned int i;
    
    switch(rotation){
        case 0:
        {
            PairCopulaVfun(family, theta, U, V, u, n);
            break;
        }
        case 90:
        {
            PairCopulaHfun(family, theta, V, U, u, n);
            break;
        }
        case 180:
        {
            PairCopulaVfun(family, theta, U, V, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
        case 270:
        {
            PairCopulaHfun(family, theta, V, U, u, n);
            for (i=0;i<n;i++) u[i]=1-u[i];
            break;
        }
    }
    
    return;
    
}

void PairCopulaVfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    switch(family){
        case 2:
        {
            unsigned int i;
            // AsymFGM
            double UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = *theta*UU*(3*UU-2)*pow(VV,3)+pow(VV,2)**theta*(1-2*UU)+1;
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 16:
        {
            // Tawn1
            PairCopulaHfun(17, theta, V, U, u, n);
            break;
        }
        case 17:
        {
            // Tawn2
            PairCopulaHfun(16, theta, V, U, u, n);
            break;
        }
        case 18:
        {
            unsigned int i;
            // Tawn
            double h1,h2,h3,h4,h5,h6,h7,hv,hu,hhv,UU,VV,uu;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            h5 = 1-theta[2];
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,theta[1]);
                hv = pow(VV,theta[2]);
                hhv = pow(VV,h5);
                h6 = pow(-log(hu),theta[0])+pow(-log(hv),theta[0]);
                h7 = exp(-pow(h6,h3));
                
                uu = h4/hu*hhv*h7 + theta[1]*hhv*pow(-log(hu),h1)/hu*(pow(h6,h2))*h7;
                *(u+i) = CheckBounds(uu);
            }
            break;
            
        }
        default:
        {
            // symmetric copulas
            PairCopulaHfun(family, theta, V, U, u, n);
            break;
        }
    }
    return;
}

void PairCopulaHfun(int family, const double *theta, double *U, double *V, double *u, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 0:
        {
            // Indep
            for (i=0;i<n;i++)
            {
                *(u+i) = U[i];
            }
            break;
        }
        case 1:
        {
            // AMH
            double h1,UU,VV,uu;
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h1 = *theta * (1-UU);
                
                uu = UU*(1-h1)/pow(1-h1*(1-VV),2);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 2:
        {
            // AsymFGM
            double UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = UU*(1+*theta*pow(1-UU,2)*(2-3*VV)*VV);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 3:
        {
            // BB1
            double h1,h3,h4,h5,h6,h7,h8,hu,hv,huv,hhuv,UU,VV,uu;
            
            h1 = -theta[0];
            h3 = 1/h1-1;
            h4 = 1/ theta[1];
            h5 = h4-1;
            h6 = theta[1] -1;
            h7 = h1-1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h8 = pow(VV,h1)-1;
                hu = pow(pow(UU,h1)-1,theta[1]);
                hv = pow(h8,theta[1]);
                huv = hu + hv;
                hhuv = 1+pow(huv,h4);
                
                uu = pow(h8,h6)*pow(VV,h7)*pow(hhuv,h3)*pow(huv,h5);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 4:
        {
            // BB6
    double h1,h2,h3,h4,h5,hu,hv,huv,hhuv,hhhuv,UU,VV,uu;
    
        h1 = theta[0]-1;
        h2 = theta[1]-1;
        h3 = 1/ theta[1];
        h4 = h3-1;
        h5 = 1/ theta[0]-1;
        
        for (i=0;i<n;i++)
        {
            UU = CheckBounds(U[i]);
            VV = CheckBounds(V[i]);
            hu = -log(1-pow(1-UU,theta[0]));
            hv = -log(1-pow(1-VV,theta[0]));
            huv = pow(hu,theta[1])+pow(hv,theta[1]);
            hhuv = exp(-pow(huv,h3));
            hhhuv = 1-hhuv;
            
            uu = pow(1-VV,h1)/(1-pow(1-VV,theta[0]))*pow(hv,h2)*pow(huv,h4)*hhuv*pow(hhhuv,h5);
            *(u+i) = CheckBounds(uu);
        }
            break;
        }
        case 5:
        {
            // BB7
            double h1,h2,h3,h4,h5,h6,hu,hv,huv,UU,VV,uu;
            
            h1 = theta[0]-1;
            h2 = 1/ theta[0]-1;
            h3 = -theta[1];
            h4 = 1/h3;
            h5 = h4 -1;
            h6 = h3-1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = 1-UU;
                hv = 1-VV;
                huv = pow(1-pow(hu,theta[0]),h3)+pow(1-pow(hv,theta[0]),h3)-1;
                
                uu = pow(1-pow(huv,h4),h2)*pow(huv,h5)*pow(1-pow(hv,theta[0]),h6)*pow(hv,h1);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 6:
        {
            // BB8
            double h1,h2,h3,h4,hu,hv,hu2,UU,VV,uu;
            
            h1 = 1-theta[1];
            h2 = 1/(1-pow(h1,theta[0]));
            h3 = 1/ theta[0]-1;
            h4 = theta[0]-1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = 1-pow(1-theta[1]*UU,theta[0]);
                hv = 1-pow(1-theta[1]*VV,theta[0]);
                hu2 = h2*hu;
                
                uu = pow(1-hu2*hv,h3)*hu2*pow(1-theta[1]*VV,h4);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 7:
        {
            //Clayton
            double h1,UU,VV,uu;
            
            h1 = -(1+*theta)/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = pow(pow(VV,*theta)*(pow(UU,-*theta)-1)+1,h1);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 8:
        {
            // FGM
            double UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = UU*(1+*theta*(1-UU)*(1-2*VV));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 9:
        {
            // Frank
            double h1,h2,h3,h4,UU,VV,uu;
            
            h1 = exp(-*theta);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h2 = pow(h1,UU);
                h3 = pow(h1,VV);
                h4 = h2 * h3;
                
                uu = (h4-h3)/(h4-h2-h3+h1);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 10:
        {
            // Gaussian
            double x, y, h1,UU,VV,uu;
            
            h1 = sqrt(1-pow(theta[0],2));
            
            boost::math::normal dist(0,1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                x = boost::math::quantile(dist, UU);
                y = boost::math::quantile(dist, VV);
                uu = boost::math::cdf(dist, (x-theta[0]*y)/h1);
                
                //x = gsl_cdf_ugaussian_Pinv(UU);
                //y = gsl_cdf_ugaussian_Pinv(VV);
                //uu = gsl_cdf_ugaussian_P((x-theta[0]*y)/h1);
                
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 11:
        {
            // Gumbel
            double h1, h2, h3, h4, h5, UU, VV, uu;
            
            h1 = *theta-1;
            h2 = (1-*theta)/ *theta;
            h3 = 1/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h4 = -log(VV);
                h5 = pow(-log(UU),*theta)+pow(h4,*theta);
                
                uu = pow(h4,h1)/VV*(pow(h5,h2))*exp(-pow(h5,h3));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 12:
        {
            // IteratedFGM
            double UU,VV,uu;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = UU*(1-theta[0]*(1-UU)*(3*theta[1]*UU*pow(VV,2)+2*VV*(1-theta[1]*UU)-1));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 13:
        {
            // Joe
            double h1,h2,h3,h4,h5,h6,UU,VV,uu;
            
            h1 = (1-*theta)/ *theta;
            h2 = *theta-1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h3 = pow(1-UU,*theta);
                h4 = 1-VV;
                h5 = pow(h4,*theta);
                h6 = h3 + h5 - h3 * h5;
                
                uu = pow(h6,h1)*pow(h4,h2)*(1-h3);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 14:
        {
            // PartialFrank
            double h1,h2,h3,hu,pUV,sUV,sMp,sMp1,UU,VV,uu;
            h1 = expm1(-*theta);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                pUV = UU*VV;
                sUV = UU+VV;
                sMp = sUV-pUV;
                sMp1 = sMp-1;
                h2 = *theta+log(1-h1*sMp1);
                h3 = sMp**theta;
                hu = 1-UU;
                
                uu = (UU*h2)/h3-(pUV*h2*hu/pow(sMp,2))/ *theta+pUV*h1*hu/((h1*sMp1-1.0)*h3);
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 15:
        {
            // Plackett
            double h1,h2,h3,h4,UU,VV,uu;
            
            h1 = *theta-1;
            h2 = 2*h1;
            h3 = 4**theta*h1;
            h4 = *theta + 1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                uu = (h1-(h2*(1+VV*h1-UU*h4))/(2*sqrt(pow(1+h1*(UU+VV),2)-h3*UU*VV)))/h2;
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 16:
        {
            // Tawn1
            double h1,h2,h3,h4,h5,hu,UU,VV,uu;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,theta[1]);
                h5 = pow(-log(hu),theta[0])+pow(-log(VV),theta[0]);
                
                uu = pow(UU,h4)*pow(-log(VV),h1)/VV*(pow(h5,h2))*exp(-pow(h5,h3));
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 17:
        {
            // Tawn2
            double h1,h2,h3,h4,h5,h6,hv,UU,VV,uu;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hv = pow(VV,theta[1]);
                h5 = pow(-log(UU),theta[0])+pow(-log(hv),theta[0]);
                h6 = exp(-pow(h5,h3));
                
                uu = h4/hv*h6+theta[1]*pow(-log(hv),h1)/hv*(pow(h5,h2))*h6;
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        case 18:
        {
            // Tawn
            double h1,h2,h3,h4,h5,h6,h7,hv,hu,hhu,UU,VV,uu;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            h5 = 1-theta[2];
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,theta[1]);
                hhu = pow(UU,h4);
                hv = pow(VV,theta[2]);
                h6 = pow(-log(hu),theta[0])+pow(-log(hv),theta[0]);
                h7 = exp(-pow(h6,h3));
                
                uu = h5/hv*hhu*h7 + theta[2]*hhu*pow(-log(hv),h1)/hv*(pow(h6,h2))*h7;
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
                
                x = boost::math::quantile(dist1, UU);
                y = boost::math::quantile(dist1, VV);
                uu = boost::math::cdf(dist2, (x-theta[0]*y)/sqrt((theta[1]+pow(y,2))*h1/nu1));
                
                //x = gsl_cdf_tdist_Pinv (UU,theta[1]);
                //y = gsl_cdf_tdist_Pinv (VV,theta[1]);
                //uu = gsl_cdf_tdist_P ((x-theta[0]*y)/sqrt((theta[1]+pow(y,2))*h1/nu1),nu1);
                
                *(u+i) = CheckBounds(uu);
            }
            break;
        }
        
    }
    
    return;
}
