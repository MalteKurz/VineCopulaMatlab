#include "VineCPP_header.hpp"

void PairCopulaPDF_Rotated_Obs(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n)
{     
    switch(rotation){
        case 0: case 180:
        {
            PairCopulaPDF(family, theta, U, V, p, n);
            break;
        }
        case 90: case 270:
        {
            PairCopulaPDF(family, theta, V, U, p, n);
            break;
        }
    }
    
    return;
}

void PairCopulaPDF(int family, int rotation, const double *theta, double *U, double *V, double *p, unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    switch(rotation){
        case 0: case 180:
        {
            PairCopulaPDF(family, theta, U, V, p, n);
            break;
        }
        case 90: case 270:
        {
            PairCopulaPDF(family, theta, V, U, p, n);
            break;
        }
    }
    
    return;
}

void PairCopulaPDF(int family, const double *theta, double *U, double *V, double *p, unsigned int n)
{
    unsigned int i;
    
    switch(family){
        case 0:
        {
            // Indep
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            break;
        }
        case 1:
        {
            // AMH
            double h1,h2,h3,h4,UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h1 = *theta * (1-UU) * (1-VV);
                h2 = (1+UU) * (1+VV);
                h3 = 1-h1;
                h4 = 1+*theta*(h2-3)+*theta*h1;
                
                p[i] = h4/pow(h3,3);
            }
            break;
        }
        case 2:
        {
            // AsymFGM
            double UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                p[i] = 1-*theta*VV*(3*VV-2)*(3*pow(UU,2)-4*UU+1);
            }
            break;
        }
        case 3:
        {
            // BB1
            if (theta[0] == 0 || theta[1] == 1)
            {
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            }
            else
            {
            //double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,hu,hv,huv,hhuv,UU,VV;
            double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10;
            
            h1 = -theta[0];
            h2 = 1+theta[0];
            h3 = 1/h1-2;
            h4 = 1/ theta[1];
            h5 = 2*h4-2;
            h6 = theta[1] -1;
            h7 = h1-1;
            h8 = h3 + 1;
            h9 = h4-2;
            h10 = theta[0]*h6;
            
            for (i=0;i<n;i++)
            {
                double h11,h12,hu,hv,huv,hhuv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h11 = pow(UU,h1)-1;
                h12 = pow(VV,h1)-1;
                hu = pow(h11,theta[1]);
                hv = pow(h12,theta[1]);
                huv = hu + hv;
                hhuv = 1+pow(huv,h4);
                
                p[i] = pow(h11,h6)*pow(h12,h6)*pow(UU,h7)*pow(VV,h7)*(h2*pow(hhuv,h3)*pow(huv,h5)+(h10*pow(hhuv,h8)*pow(huv,h9)));
            }
            }
            break;
        }
        case 4:
        {
            // BB6
            //double h1,h2,h3,h4,h5,h6,u1,v1,ut,vt,hu,hv,huv,hhuv,ehhuv,hhhuv,UU,VV;
            double h1,h2,h3,h4,h5,h6;
            
            h1 = theta[0]-1;
            h2 = theta[1]-1;
            h3 = 1/ theta[1];
            h4 = h3-2;
            h5 = 1/ theta[0]-1;
            h6 = theta[0]*h2;
            
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
                hhuv = pow(huv,h3);
                ehhuv = exp(-hhuv);
                hhhuv = 1-ehhuv;
                
                p[i] = pow(u1,h1)*pow(v1,h1)/ut/vt*pow(hu,h2)*pow(hv,h2)*pow(huv,h4)*ehhuv*pow(hhhuv,h5)*(h1*hhuv*ehhuv/hhhuv+theta[0]*hhuv+h6);
            }
            break;
        }
        case 5:
        {
            // BB7
            //double h1,h2,h3,h4,h5,h6,h7,hu,hv,huv,UU,VV;
            double h1,h2,h3,h4,h5,h6,h7;
            
            h1 = theta[0]-1;
            h2 = theta[0]*(theta[1]+1);
            h3 = 1/ theta[0]-2;
            h4 = -theta[1];
            h5 = -2/ theta[1] -2;
            h6 = 1/ theta[1];
            h7 = h4-1;
            
            for (i=0;i<n;i++)
            {
                double hu,hv,huv,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = 1-UU;
                hv = 1-VV;
                huv = pow(1-pow(hu,theta[0]),h4)+pow(1-pow(hv,theta[0]),h4)-1;
                
                p[i] = pow(1-pow(hu,theta[0]),h7)*pow(hu,h1)*pow(1-pow(hv,theta[0]),h7)*pow(hv,h1)*pow(1-pow(huv,-h6),h3)*pow(huv,h5)*(h1+h2*(pow(huv,h6)-1));
            }
            break;
        }
        case 6:
        {
            // BB8
            if (theta[0] == 1 || theta[1] == 0)
            {
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            }
            else
            {
                //double h1,h2,h3,h4,h5,h6,hu,hv,huv,UU,VV;
                double h1,h2,h3,h4,h5;
                
                h1 = 1-theta[1];
                h2 = pow(h1,theta[0]);
                h3 = 1/(1-h2);
                h4 = 1/ theta[0];
                h5 = -theta[0]*h2+theta[0]-1;
                
                for (i=0;i<n;i++)
                {
                    double h6,hu,hv,huv,UU,VV;
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    hu = pow(1-theta[1]*UU,theta[0]);
                    hv = pow(1-theta[1]*VV,theta[0]);
                    huv = hu*hv-hu-hv;
                    
                    p[i] = theta[1]*hu*hv*pow(1-h3*(1-hu)*(1-hv),h4)*(h5-huv)/((theta[1]*UU-1)*(theta[1]*VV-1))/(pow(h2+huv,2));
                }
            }
            break;
        }
        case 7:
        {
            //Clayton
            if (*theta == 0)
            {
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            }
            else
            {
                double h1,h2,h3,h4,UU,VV;
                
                h1 = (1+2**theta)/ *theta;
                h2 = 1+*theta;
                h3 = -*theta;
                
                for (i=0;i<n;i++)
                {
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    
                    h4 = pow(UU,h3)+pow(VV,h3)-1;
                    
                    p[i] = h2*pow(UU,-h2)*pow(VV,-h2)*pow(h4,-h1);
                }
            }
            break;
        }
        case 8:
        {
            // FGM
            double UU,VV;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                p[i] = 1+*theta*(1-2*UU)*(1-2*VV);
            }
            break;
        }
        case 9:
        {
            // Frank
            if (*theta == 0)
            {
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            }
            else
            {
                double h1,h2,h3,h4,UU,VV;
                
                h1 = -*theta;
                h2 = expm1(h1);
                h3 = h1*h2;
                
                for (i=0;i<n;i++)
                {
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    h4 = expm1(h1*UU)*expm1(h1*VV);
                    
                    p[i] = h3*exp(h1*(UU+VV))/pow(h2+h4,2);
                }
            }
            break;
        }
        case 10:
        {
            // Gaussian
            double x,y,rho2,h1,h2,h3,UU,VV;
            
            rho2 = pow(*theta,2);
            h1 = 1-rho2;
            h2 = rho2/(2*h1);
            h3 = *theta /h1;
            
            boost::math::normal dist(0,1);
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                x = boost::math::quantile(dist, UU);
                y = boost::math::quantile(dist, VV);
                
                //x = gsl_cdf_ugaussian_Pinv(UU);
                //y = gsl_cdf_ugaussian_Pinv(VV);

                p[i] = exp(h3*x*y-h2*(pow(x,2)+pow(y,2)))/sqrt(h1);
            }
            break;
        }
        case 11:
        {
            // Gumbel
            double h1,h2,h3,h4,h5,h6,h7,UU,VV;
            
            h1 = *theta-1;
            h2 = (1-2**theta)/ *theta;
            h3 = 1/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                h4 = -log(UU);
                h5 = -log(VV);
                h6 = pow(h4,*theta) + pow(h5,*theta);
                h7 = pow(h6,h3);
                
                p[i] = exp(-h7+h4+h5)*pow(h4,h1)*pow(h5,h1)*pow(h6,h2)*(h1+h7);
            }
            break;
        }
        case 12:
        {
            // IteratedFGM
            double h1,h2,h3,h4,h5,h6,pUV,sUV,UU,VV;
            h1 = 1+theta[0];
            h2 = -2*theta[0];
            h3 = theta[0]*theta[1];
            h4 = 4*(h3+theta[0]);
            h5 = 9*h3;
            h6 = 6*h3;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                sUV = UU + VV;
                pUV = UU*VV;
                
                p[i] = h1+h2*sUV+pUV*(h4+pUV*h5-sUV*h6);
            }
            break;
        }
        case 13:
        {
            // Joe
            double h1,h2,h3,h4,h5,h6,u1,v1,h31,h41,UU,VV;
            
            h1 = (1-2**theta)/ *theta;
            h2 = *theta-1;
            h6 = (1-*theta)/ *theta;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                u1 = 1-UU;
                h3 = pow(u1,*theta);
                h31 = pow(u1,h2);
                v1 = 1-VV;
                h4 = pow(v1,*theta);
                h41 = pow(v1,h2);
                h5 = h3 + h4 - h3 * h4;
                
                p[i] = h2*pow(h5,h1)*h31*(1-h4)*h41*(1-h3)+pow(h5,h6)**theta*h31*h41;
            }
            break;
        }
        case 14:
        {
            // PartialFrank
            if (*theta == 0)
            {
            for (i=0;i<n;i++)
            {
                p[i] = 1;
            }
            }
            else
            {
                double h1,h2,h3,h4,h5,h6,h4_2,UU,VV;
                
                h1 = expm1(-*theta);
                
                for (i=0;i<n;i++)
                {
                    UU = CheckBounds(U[i]);
                    VV = CheckBounds(V[i]);
                    h2 = UU*VV;
                    h3 = 2*h2;
                    h4 = UU+VV-h2;
                    h5 = h4-1;
                    h4_2 = pow(h4,2);
                    
                    h6 = h1*h5-1;
                    
                    p[i] = (h1*(h4_2-h3)*h4/h6+h5*h2*h4_2*pow(h1,2)/pow(h6,2)+h3*(*theta+log(-h1*h5+1)))/(*theta*pow(h4,3));
                }
            }
            break;
        }
        case 15:
        {
            // Plackett
            double h1,h2,huv,hhuv,UU,VV;
            
            h1 = *theta-1;
            h2 = 4**theta*h1;
            
            for (i=0;i<n;i++)
            {
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                huv = UU+VV;
                hhuv = UU*VV;
                
                p[i] = *theta*(h1*(huv-2*hhuv)+1)/pow(pow(1+h1*huv,2)-h2*hhuv,1.5);
            }
            break;
        }
        case 16:
        {
            // Tawn1
            //double h1,h2,h3,h4,h5,hu,h2g,h4g,h5g,h6g,h7g,UU,VV;
            double h1,h2,h3,h4,h5,h2g;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            h5 = -theta[1];
            h2g = (1-2*theta[0])/ theta[0];
            
            for (i=0;i<n;i++)
            {
                double hu,h4g,h5g,h6g,h7g,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,theta[1]);
                h4g = -log(hu);
                h5g = -log(VV);
                h6g = pow(h4g,theta[0]) + pow(h5g,theta[0]);
                h7g = pow(h6g,h3);
                
                p[i] = h4*pow(UU,h5)*pow(h5g,h1)/VV*(pow(h6g,h2))*exp(-pow(h6g,h3))+theta[1]*(exp(-h7g)/hu/VV*pow(h4g,h1)*pow(h5g,h1)*pow(h6g,h2g)*(h1+h7g));
            }
            break;
        }
        case 17:
        {
            // Tawn2
            //double h1,h2,h3,h4,h5,hv,h2g,h4g,h5g,h6g,h7g,UU,VV;
            double h1,h2,h3,h4,h5,h2g;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h4 = 1-theta[1];
            h5 = -theta[1];
            h2g = (1-2*theta[0])/ theta[0];
            
            for (i=0;i<n;i++)
            {
                double hv,h4g,h5g,h6g,h7g,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hv = pow(VV,theta[1]);
                h4g = -log(hv);
                h5g = -log(UU);
                h6g = pow(h4g,theta[0]) + pow(h5g,theta[0]);
                h7g = pow(h6g,h3);
                
                p[i] = h4*pow(VV,h5)*pow(h5g,h1)/UU*(pow(h6g,h2))*exp(-pow(h6g,h3))+theta[1]*(exp(-h7g)/hv/UU*pow(h4g,h1)*pow(h5g,h1)*pow(h6g,h2g)*(h1+h7g));
            }
            break;
        }
        case 18:
        {
            // Tawn
            //double h1,h2,h3,h4,hu,hv,h2g,h41,h42,h412,h41d,h42d,h4g,h4g1,h5g,h5g1,h6g,h7g,h8g,h9g,h10g,UU,VV;
            double h1,h2,h3,h4,h2g,h41,h42,h412,h41d,h42d;
            
            h1 = theta[0]-1;
            h2 = (1-theta[0])/ theta[0];
            h3 = 1/ theta[0];
            h41 = (1-theta[1]);
            h42 = (1-theta[2]);
            h412 = h41*h42;
            h41d = h41*theta[2];
            h42d = h42*theta[1];
            h2g = (1-2*theta[0])/ theta[0];
            
            for (i=0;i<n;i++)
            {
                double hu,hv,h4g,h4g1,h5g,h5g1,h6g,h7g,h8g,h9g,h10g,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                hu = pow(UU,theta[1]);
                hv = pow(VV,theta[2]);
                h4g = -log(hu);
                h5g = -log(hv);
                h4g1 = pow(h4g,h1);
                h5g1 = pow(h5g,h1);
                h6g = pow(h4g,theta[0]) + pow(h5g,theta[0]);
                h7g = pow(h6g,h3);
                h8g = exp(-h7g);
                h9g = h8g/hu/hv;
                h10g = pow(h6g,h2);
                
                p[i]= h9g*(theta[1]*theta[2]*(h4g1*h5g1*pow(h6g,h2g)*(h1+h7g))
                +h41d*h5g1*h10g
                        +h42d*h4g1*h10g
                        +h412);
            }
            break;
        }
        case 19:
        {
            // t
            //double x,y,x2,y2,rho2,h1,h2,h3,h4,h5,h6,UU,VV;
            double rho2,h1,h2,h3,h4,h5,h6;
            
            rho2 = pow(theta[0],2);
            h1 = 1-rho2;
            h2 = theta[1]/2;
            h3 = h2 + 0.5;
            h4 = h2 + 1;
            h5 = 1/ theta[1];
            h6 = h5/h1;
            
            boost::math::students_t dist(theta[1]);
            
            for (i=0;i<n;i++)
            {
                double x,y,x2,y2,UU,VV;
                UU = CheckBounds(U[i]);
                VV = CheckBounds(V[i]);
                
                x = boost::math::quantile(dist, UU);
                y = boost::math::quantile(dist, VV);
                
                //x = gsl_cdf_tdist_Pinv (UU, theta[1]);
                //y = gsl_cdf_tdist_Pinv (VV, theta[1]);
                
                x2 = pow(x,2);
                y2 = pow(y,2);
                
                p[i] = tgamma(h4)*tgamma(h2)/sqrt(h1)/pow(tgamma(h3),2)*pow(1+h5*x2,h3)*pow(1+h5*y2,h3)/pow(1+h6*(x2+y2-2*theta[0]*x*y),h4);
            }
            break;
        }
        
    }
    
    return;
}
