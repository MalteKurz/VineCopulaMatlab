#include "VineCPP_header.hpp"

typedef struct {
    int family;
    double *U;
    double *V;
    unsigned int n;
} my_data;

double ObjectiveFunction(unsigned n, const double *x, double *grad, void *data)
{
    my_data *data_for_func = (my_data *) data;
    double CLL = PairCopulaNegLL(data_for_func->family,x,data_for_func->U,data_for_func->V,data_for_func->n);
    return CLL;
}

void PairCopulaFit_Rotated_Obs(double *theta,int family, int rotation, double *U,double *V,unsigned int n)
{     
    switch(rotation){
        case 0: case 180:
        {
            return PairCopulaFit(theta, family, U, V, n);
            break;
        }
        case 90: case 270:
        {
            return PairCopulaFit(theta, family, V, U, n);
            break;
        }
    }
}

void PairCopulaFit(double *theta,int family, int rotation, double *U,double *V,unsigned int n)
{     
    if(rotation>0)
    {
        Rotate_Obs(U,V,rotation,n);
    }
    
    switch(rotation){
        case 0: case 180:
        {
            return PairCopulaFit(theta, family, U, V, n);
            break;
        }
        case 90: case 270:
        {
            return PairCopulaFit(theta, family, V, U, n);
            break;
        }
    }
}

void PairCopulaFit(double *theta,int family, double *U,double *V,unsigned int n)
{
    if (family == 0)
    {
        return;
    }
    else
    {
        std::vector<double> bounds(120);
        LoadBounds(&bounds[0]);
        nlopt::opt opt;
        std::vector<double> lb(1);
        std::vector<double> ub(1);
        std::vector<double> x(1);
        
        
        // Depending on the family associate outputs and upper and lower bounds
        switch(family){
            case 18:
            {
                opt = nlopt::opt(nlopt::LN_BOBYQA,3);
                lb.resize(3);
                ub.resize(3);
                x.resize(3);
                break;
            }
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                opt = nlopt::opt(nlopt::LN_BOBYQA,2);
                lb.resize(2);
                ub.resize(2);
                x.resize(2);
                break;
            }
            default:
            {
                opt = nlopt::opt(nlopt::LN_BOBYQA,1);
            }
        }
        
        
        // Getting the lower and upper bounds
        switch(family){
            case 1:
            {
                // AMH
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 0;
                break;
            }
            case 2:
            {
                // AsymFGM
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 0.5;
                break;
            }
            case 3:
            {
                // BB1
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 0.5;
                x[1] = 2;
                break;
            }
            case 4:
            {
                // BB6
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 2;
                x[1] = 2;
                break;
            }
            case 5:
            {
                // BB7
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 2;
                x[1] = 0.25;
                break;
            }
            case 6:
            {
                // BB8
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 2;
                x[1] = 0.25;
                break;
            }
            case 7:
            {
                //Clayton
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                break;
            }
            case 8:
            {
                // FGM
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 0;
                break;
            }
            case 9:
            {
                // Frank
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                break;
            }
            case 10:
            {
                // Gaussian
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 0;
                break;
            }
            case 11:
            {
                // Gumbel
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                break;
            }
            case 12:
            {
                // IteratedFGM
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 0;
                x[1] = 0;
                break;
            }
            case 13:
            {
                // Joe
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                break;
            }
            case 14:
            {
                // PartialFrank
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                break;
            }
            case 15:
            {
                // Plackett
                lb[0] = bounds[family*6];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                break;
            }
            case 16:
            {
                // Tawn1
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                x[1] = 0.5;
                break;
            }
            case 17:
            {
                // Tawn2
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 5;
                x[1] = 0.5;
                break;
            }
            case 18:
            {
                // Tawn
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                lb[2] = bounds[family*6+4];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                ub[2] = bounds[family*6+5];
                opt.set_upper_bounds(ub);
                x[0] = 2;
                x[1] = 0.5;
                x[2] = 0.5;
                break;
            }
            case 19:
            {
                // t
                lb[0] = bounds[family*6];
                lb[1] = bounds[family*6+2];
                opt.set_lower_bounds(lb);
                ub[0] = bounds[family*6+1];
                ub[1] = bounds[family*6+3];
                opt.set_upper_bounds(ub);
                x[0] = 0.5;
                x[1] = 5;
                break;
            }
            
        }
        
        my_data data = {family,U,V,n};
        
        opt.set_min_objective(ObjectiveFunction,&data);
        
        opt.set_xtol_rel(1e-5);
        
        double minf;
        try
        {
            //nlopt::result result = opt.optimize(x, minf);
            opt.optimize(x, minf);
        }
        catch(nlopt::roundoff_limited) {
            //mexWarnMsgTxt("Halted because roundoff errors limited progress.");
        }

        
        //mexPrintf("result code = %i \n", opt.last_optimize_result());
        
        void nlopt_destroy(nlopt_opt opt);
        
        switch(family){
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                theta[0] = x[0];
                theta[1] = x[1];
                break;
            }
            case 18:
            {
                theta[0] = x[0];
                theta[1] = x[1];
                theta[2] = x[2];
                break;
            }
            default:
            {
                theta[0] = x[0];
                break;
            }
        }
    }
    
    return;
    
}
