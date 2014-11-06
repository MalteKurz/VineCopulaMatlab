#include "VineCPP_header.hpp"

int UpdateStructure(std::vector<int>& structure, int StructuringRule, double *U, unsigned int d, unsigned int d1, unsigned int n)
{
    unsigned int i,j;
    
    std::vector<double> tau(d1*d1);
    SD_Kendall_Tau_Matrix(&tau[0],U,d1,n);
    
    std::vector<double> tau1(d1);
    double Tau1;
    
    for (i=0;i<d1;i++)
    {
        for (j=0;j<d1;j++)
        {
            Tau1 = std::abs(tau[i*d1+j]);
            tau1[i] += Tau1;
        }
    }
    
    std::vector<double>::iterator Maximum;
    
    Maximum = std::max_element(tau1.begin(), tau1.end());
    int I = std::distance(tau1.begin(), Maximum);
    
    int idxI = structure[I+d-d1];
    
    structure.erase (structure.begin()+I+d-d1);
    
    structure.insert (structure.begin()+d-d1,idxI);
    switch (StructuringRule)
    {
        case 1:
        {
            tau1.resize(d1-1);
            for (j=0;j<I;j++)
            {
                Tau1 = std::abs(tau[I*d1+j]);
                tau1[j] = Tau1;
            }
            for (j=I+1;j<d1;j++)
            {
                Tau1 = std::abs(tau[I*d1+j]);
                tau1[j-1] = Tau1;
            }
            
            for (j=0;j<d1-2;j++)
            {
                Maximum = std::max_element(tau1.begin(), tau1.end());
                int I = std::distance(tau1.begin(), Maximum);
                
                tau1.erase(tau1.begin()+I);
                
                int idxI = structure[I+1+j];
                
                structure.erase (structure.begin()+I+1+j);
                
                structure.insert (structure.begin()+1+j,idxI);
            }
            break;
        }
        case 2:
        {
            tau1.resize(d1-1);
            for (j=0;j<I;j++)
            {
                Tau1 = std::abs(tau[I*d1+j]);
                tau1[j] = Tau1;
            }
            for (j=I+1;j<d1;j++)
            {
                Tau1 = std::abs(tau[I*d1+j]);
                tau1[j-1] = Tau1;
            }
            
            for (j=0;j<d1-2;j++)
            {
                Maximum = std::min_element(tau1.begin(), tau1.end()); // Note that Maximum is the minimum in this case!
                int I = std::distance(tau1.begin(), Maximum);
                
                tau1.erase(tau1.begin()+I);
                
                int idxI = structure[I+1+j];
                
                structure.erase (structure.begin()+I+1+j);
                
                structure.insert (structure.begin()+1+j,idxI);
            }
            break;
        }
        case 3:
        {
            for (j=0;j<d1-2;j++)
            {
                tau1.resize(d1-1-j);
                for (i=0;i<d1-1-j;i++)
                {
                    Tau1 = std::abs(tau[(I+j)*d1+structure[i+1+j]]);
                    tau1[i] = Tau1;
                }
                Maximum = std::min_element(tau1.begin(), tau1.end()); // Note that Maximum is the minimum in this case!
                int I = std::distance(tau1.begin(), Maximum);
                
                int idxI = structure[I+1+j];
                
                structure.erase (structure.begin()+I+1+j);
                
                structure.insert (structure.begin()+1+j,idxI);
            }
            break;
        }
        case 0:
        {
            
        }
    }
    return I;
}

void TreeSelect(int type, int *families, std::vector<double>& thetas, int *rotations, double *U, int d, unsigned int n, double *familyset, int m)
{
    int i;
    
    std::vector<double> theta(3*(d-1));
    
    switch (type)
    {
        case 0: // C-Vine (Pair-copula selection for the whole vine)
        {
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                PairCopulaSelect(&families[i], &theta[3*i], &rotations[i], &U[0], &U[(i+1)*n], n, familyset, m);
            }
            break;
        }
        case 1: // D-Vine (Pair-copula selection for the first tree)
        {
            #pragma omp parallel for private(i)
            for (i=0;i<d-1;i++)
            {
                PairCopulaSelect(&families[i], &theta[3*i], &rotations[i], &U[i*n], &U[(i+1)*n], n, familyset, m);
            }
        }
    }
    
    for (i=0;i<d-1;i++)
    {
        switch(families[i]){
            case 0:
            {
                break;
            }
            case 18:
            {
                thetas.reserve(thetas.size() + 3);
                thetas.insert(thetas.end(), &theta[3*i], &theta[3*i+3]);
                break;
            }
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                thetas.reserve(thetas.size() + 2);
                thetas.insert(thetas.end(), &theta[3*i], &theta[3*i+2]);
                break;
            }
            default:
            {
                thetas.push_back(theta[3*i]);
            }
        }
    }
}

// Tree selection for higher (>2) order trees of D-Vine copulas
void TreeSelect(int type, int *families, std::vector<double>& thetas, int *rotations, double *H, double *V, int d, unsigned int n, double *familyset, int m)
{
    int i;
    
    std::vector<double> theta(3*(d-1));
    
    
    #pragma omp parallel for private(i)
    for (i=0;i<d-1;i++)
    {
        PairCopulaSelect(&families[i], &theta[3*i], &rotations[i], &H[i*n], &V[i*n], n, familyset, m);
    }
    
    for (i=0;i<d-1;i++)
    {
        switch(families[i]){
            case 0:
            {
                break;
            }
            case 18:
            {
                thetas.reserve(thetas.size() + 3);
                thetas.insert(thetas.end(), &theta[3*i], &theta[3*i+3]);
                break;
            }
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                thetas.reserve(thetas.size() + 2);
                thetas.insert(thetas.end(), &theta[3*i], &theta[3*i+2]);
                break;
            }
            default:
            {
                thetas.push_back(theta[3*i]);
            }
        }
    }
}

void GetPseudoObs(int *families, std::vector<double>& thetas, int *rotations, double *V, int d, unsigned int n)
{
    std::vector<int> NumbParams(d-1);
    std::vector<double> U1(n),V1(n);
    
    unsigned int i;
    int J=0;
    for (i=d-1;i>=1;i--)
    {
        switch(families[i-1]){
            case 0:
            {
                break;
            }
            case 18:
            {
                J += 3;
                break;
            }
            case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
            {
                J += 2;
                break;
            }
            default:
            {
                J += 1;
            }
        }
        NumbParams[i-1] = J;
    }
    
    int TotalNumbParams = thetas.size();
    
    for (i=d-1;i>=1;i--)
    {
        if (families[i-1] != 0)
        {
            if (rotations[i-1]>0)
            {
                Rotate_Obs(&V[0], &V[i*n],&U1[0],&V1[0],rotations[i-1],n);
                PairCopulaVfun_Rotated_Obs(families[i-1], rotations[i-1], &thetas[TotalNumbParams-NumbParams[i-1]], &U1[0], &V1[0] , &V[i*n], n);
            }
            else
            {
                PairCopulaVfun(families[i-1], &thetas[TotalNumbParams-NumbParams[i-1]], &V[0], &V[i*n], &V[i*n], n);
            }
        }
    }
}

void VineCopulaStructureSelect(int type, double *Structure, double *Families, double *Rotations, std::vector<double>& Thetas, double *U, unsigned int d, unsigned int n, int StructuringRule, double *familyset, int m)
{
    std::vector<int> families((d-1)*d/2);
    std::vector<int> rotations((d-1)*d/2);
    std::vector<double> thetas;
    
    int i, j, k, I;
    
    switch (type)
    {
        case 0: // C-Vine
        {
            std::vector<int> structure(d);
            
            for (i=0;i<d;i++)
            {
                structure[i] = i;
            }
            
            switch (StructuringRule)
            {
                case 0: // Maximizing the sum of absolute Kendall's tau in every tree
                {
                    I = UpdateStructure(structure, 0, U, d, d, n);
                    
                    std::vector<double> V(n*d);
                    
                    for (i=0;i<n;i++)
                    {
                        V[i] = U[I*n+i];
                    }
                    
                    for (j=0;j<I;j++)
                    {
                        for (i=0;i<n;i++)
                        {
                            V[(j+1)*n+i] = U[j*n+i];
                        }
                    }
                    
                    for (j=I+1;j<d;j++)
                    {
                        for (i=0;i<n;i++)
                        {
                            V[j*n+i] = U[j*n+i];
                        }
                    }
                    
                    TreeSelect(type, &families[0], thetas, &rotations[0], &V[0], d, n, familyset, m);
                    
                    GetPseudoObs(&families[0], thetas, &rotations[0], &V[0], d, n);
                    
                    for (k=1;k<d-2;k++)
                    {
                        I = UpdateStructure(structure, 0, &V[n], d, d-k, n);
                        
                        for (i=0;i<n;i++)
                        {
                            V[i] = V[(I+1)*n+i];
                        }
                        
                        for (j=I+1;j<d-k;j++)
                        {
                            for (i=0;i<n;i++)
                            {
                                V[j*n+i] = V[(j+1)*n+i];
                            }
                        }
                        
                        TreeSelect(type, &families[d*k-k*(k+1)/2], thetas, &rotations[d*k-k*(k+1)/2], &V[0], d-k, n, familyset, m);
                        
                        GetPseudoObs(&families[d*k-k*(k+1)/2], thetas, &rotations[d*k-k*(k+1)/2], &V[0], d-k, n);
                    }
                    
                    TreeSelect(type, &families[d*(d-1)/2-1], thetas, &rotations[d*(d-1)/2-1], &V[n], 2, n, familyset, m);
                    break;
                }
// case 2: Choose the most dependent variable as the root node and list the other variables by their dependence to the root node in increasing order (Nikoloulopoulos et al. 2012, p.3665)
// case 3: Choose the most dependent variable as the root node and list the other variables sequentially by choosing the variable which is least dependent with the previously selcted one (Nikoloulopoulos et al. 2012, p.3665)
// case 1: Choose the most dependent variable as the root node and list the other variables by their dependence to the root node in decreasing order (Nikoloulopoulos et al. 2012, p.3665)
                default:
                {
                    I = UpdateStructure(structure, StructuringRule, U, d, d, n);
                    
                    std::vector<double> V(n*d);
                    
                    for (j=0;j<d;j++)
                    {
                        for (i=0;i<n;i++)
                        {
                            V[j*n+i] = U[structure[j]*n+i];
                        }
                    }
                    
                    TreeSelect(type, &families[0], thetas, &rotations[0], &V[0], d, n, familyset, m);
                    
                    GetPseudoObs(&families[0], thetas, &rotations[0], &V[0], d, n);
                    
                    
                    for (k=1;k<d-2;k++)
                    {
                        TreeSelect(type, &families[d*k-k*(k+1)/2], thetas, &rotations[d*k-k*(k+1)/2], &V[k*n], d-k, n, familyset, m);
                        
                        GetPseudoObs(&families[d*k-k*(k+1)/2], thetas, &rotations[d*k-k*(k+1)/2], &V[k*n], d-k, n);
                    }
                    
                    TreeSelect(type, &families[d*(d-1)/2-1], thetas, &rotations[d*(d-1)/2-1], &V[(d-2)*n], 2, n, familyset, m);
                }
            }
            
            
            Thetas.resize(thetas.size());
            
            switch (StructuringRule)
            {
                case 0:
                {
                    std::vector<int> NumbParams((d-1)*d/2+1);
                    int J=0;
                    
                    for (i=0;i<(d-1)*d/2;i++)
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
                    std::pair<int, int> R[d-1];
                    std::vector<int> Ranks(d-1);
                    
                    for (k=0;k<d-1;k++)
                    {
                        // Computing ranks
                        for (i=0;i<d-k-1;i++)
                        {
                            R[i] = std::make_pair(structure[i+1+k],i);
                        }
                        
                        std::sort(&R[0],&R[d-k-1]);
                        
                        for (i=0;i<d-k-1;i++)
                        {
                            Ranks[R[i].second] = i;
                        }
                        
                        for (i=0;i<d-k-1;i++)
                        {
                            Families[d*k-k*(k+1)/2+i] = (double) families[d*k-k*(k+1)/2+Ranks[i]];
                            Rotations[d*k-k*(k+1)/2+i] = (double) rotations[d*k-k*(k+1)/2+Ranks[i]];
                            switch(NumbParams[d*k-k*(k+1)/2+Ranks[i]+1]-NumbParams[d*k-k*(k+1)/2+Ranks[i]]){
                                case 0:
                                {
                                    break;
                                }
                                case 1:
                                {
                                    Thetas[J] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]];
                                    J++;
                                    break;
                                }
                                case 2:
                                {
                                    Thetas[J] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]];
                                    Thetas[J+1] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]+1];
                                    J = J+2;
                                    break;
                                }
                                default:
                                {
                                    Thetas[J] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]];
                                    Thetas[J+1] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]+1];
                                    Thetas[J+2] = thetas[NumbParams[d*k-k*(k+1)/2+Ranks[i]]+2];
                                    J = J+3;
                                }
                            }
                        }
                    }
                    break;
                }
                default:
                {
                    for (i=0;i<(d-1)*d/2;i++)
                    {
                        Families[i] = (double) families[i];
                        Rotations[i] = (double) rotations[i];
                    }
                    
                    for (i=0;i<(int) thetas.size();i++)
                    {
                        Thetas[i] = thetas[i];
                    }
                    
                }
            }
            
            for (i=0;i<d;i++)
            {
                Structure[i] = (double) structure[i];
            }
            break;
        }
        case 1: // D-Vine
        {
            std::vector<double> V((d-2)*n);
            std::vector<double> H((d-2)*n);
            std::vector<double> U1(n),V1(n);
            
            TreeSelect(type, &families[0], thetas, &rotations[0], &U[0], d, n, familyset, m);
            
            int J =0;
            
            for (i=0;i<d-1;i++)
            {
                if (rotations[i]>0)
                {
                    Rotate_Obs(&U[i*n],&U[(i+1)*n],&U1[0],&V1[0],rotations[i],n);
                    if (i<d-2) {
                        PairCopulaHfun_Rotated_Obs(families[i], rotations[i], &thetas[J], &U1[0], &V1[0], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun_Rotated_Obs(families[i], rotations[i], &thetas[J], &U1[0], &V1[0], &V[(i-1)*n], n);
                    }
                }
                else
                {
                    if (i<d-2) {
                        PairCopulaHfun(families[i], &thetas[J], &U[i*n], &U[(i+1)*n], &H[i*n], n);
                    }
                    if (i>0) {
                        PairCopulaVfun(families[i], &thetas[J], &U[i*n], &U[(i+1)*n], &V[(i-1)*n], n);
                    }
                }
                switch(families[i]){
                    case 0:
                    {
                        break;
                    }
                    case 18:
                    {
                        J += 3;
                        break;
                    }
                    case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
                    {
                        J += 2;
                        break;
                    }
                    default:
                    {
                        J += 1;
                    }
                }
            }

            
            for (k=1;k<d-2;k++)
            {
                TreeSelect(type, &families[d*k-k*(k+1)/2], thetas, &rotations[d*k-k*(k+1)/2], &H[0], &V[0], d-k, n, familyset, m);
                
                for (i=0;i<d-k-1;i++)
                {
                    if (rotations[d*k-k*(k+1)/2+i]>0)
                    {
                        Rotate_Obs(&H[i*n],&V[i*n],&U1[0],&V1[0],rotations[d*k-k*(k+1)/2+i],n);
                        if (i>0) {
                            PairCopulaVfun_Rotated_Obs(families[d*k-k*(k+1)/2+i], rotations[d*k-k*(k+1)/2+i], &thetas[J], &U1[0], &V1[0], &V[(i-1)*n], n);
                        }
                        if (i<d-k-2) {
                            PairCopulaHfun_Rotated_Obs(families[d*k-k*(k+1)/2+i], rotations[d*k-k*(k+1)/2+i], &thetas[J], &U1[0], &V1[0], &H[i*n], n);
                        }
                    }
                    else
                    {
                        if (i>0) {
                            PairCopulaVfun(families[d*k-k*(k+1)/2+i], &thetas[J], &H[i*n], &V[i*n], &V[(i-1)*n], n);
                        }
                        if (i<d-k-2) {
                            PairCopulaHfun(families[d*k-k*(k+1)/2+i], &thetas[J], &H[i*n], &V[i*n], &H[i*n], n);
                        }
                    }
                    switch(families[d*k-k*(k+1)/2+i]){
                        case 0:
                        {
                            break;
                        }
                        case 18:
                        {
                            J += 3;
                            break;
                        }
                        case 3: case 4: case 5: case 6: case 12: case 16: case 17: case 19:
                        {
                            J += 2;
                            break;
                        }
                        default:
                        {
                            J += 1;
                        }
                    }
                }
            }

            TreeSelect(type, &families[d*(d-1)/2-1], thetas, &rotations[d*(d-1)/2-1], &H[0], &V[0], 2, n, familyset, m);
            
            Thetas.resize(thetas.size());
            
            for (i=0;i<(d-1)*d/2;i++)
            {
                Families[i] = (double) families[i];
                Rotations[i] = (double) rotations[i];
            }
            
            for (i=0;i<(int) thetas.size();i++)
            {
                Thetas[i] = thetas[i];
            }
        }
    }
    
    return;
}
