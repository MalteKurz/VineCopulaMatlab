function [lb,ub] = PairCopulaParameterBounds(families)
%PAIRCOPULAPARAMETERBOUNDS Getting the parameter bounds for pair-copulas
% Purpose
%        The function gives back bounds for the parameter values. Those can
%        be set globally for the VineCPP toolbox by using the GUI
%        SetGlobalParameterBounds. The bounds are used for checking for
%        correctly specified pair-copulas and also as bounds for the
%        estimation of the parameters.
%
%
% Usage
%               [lb,ub] = PairCopulaParameterBounds(families)
%
%
% Inputs
%       families        = Either the copula family or a cell array or
%                         numerically coded vector of copula families.
%
%
% Outputs
%       lb              = If the input variable families is just one single
%                         copula, then it is the vector of lower bounds of
%                         the corresponding copula parameter(s). If the input
%                         variable families is a cell array of copula
%                         families, then it is a vector of lower bounds for
%                         all the copulas being specified in the families
%                         cell array.
%       ub              = If the input variable families is just one single
%                         copula, then it is the vector of upper bounds of
%                         the corresponding copula parameter(s). If the
%                         input variable families is a cell array of copula
%                         families, then it is a vector of upper bounds for
%                         all the copulas being specified in the families
%                         cell array.
%
%
%
% Author: Malte Kurz

bounds = GlobalPairCopulaParameterBounds;

if ischar(families)
    families = {families};
end

NumbPCs = length(families);

if not(isnumeric(families))
    Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
    for i=1:NumbPCs
        if sum(strcmpi(families{i},Families))
            families{i} = find(strcmp(families{i},Families))-1;
        else
            error(['The copula family ' families{i} ' is not implemented'])
        end
    end
    families = cell2mat(families);
end

NumbParams = NumbPCs + sum(families==3) + sum(families==4) + sum(families==5) + sum(families==6) + sum(families==12) + sum(families==16) + sum(families==17) + sum(families==19)...
     + 2.*sum(families==18) - sum(families==0);
ub = zeros(1,NumbParams);
lb = ub;

if length(families) == sum(families==0)
    lb = [];
    ub = [];
else
    k = 1;
    for j = 1:NumbPCs
        switch families(j)
            case {3,4,5,6,12,16,17,19}
                F = families(j)+1;
                lb(k) = bounds(F,1);
                ub(k) = bounds(F,2);
                lb(k+1) = bounds(F,3);
                ub(k+1) = bounds(F,4);
                k = k+2;
                
            case {18}
                F = families(j)+1;
                lb(k) = bounds(F,1);
                ub(k) = bounds(F,2);
                lb(k+1) = bounds(F,3);
                ub(k+1) = bounds(F,4);
                lb(k+2) = bounds(F,5);
                ub(k+2) = bounds(F,6);
                k = k+3;
                
            case {0}
                
            otherwise
                F = families(j)+1;
                lb(k) = bounds(F,1);
                ub(k) = bounds(F,2);
                k = k+1;
                
        end
    end
end

end
