function [] = CheckParameters(families,theta)
%CHECKPARAMETERS Checking for correctly specified pair-copula parameters
% Purpose:
%         The function is internally used for checking whether pair-copula
%         parameters are correctly specified. The bounds for copula
%         parameters, which are used in this function for consistency
%         checking as well as within the numerically computation of
%         maximum-likelihood estimates, can be set by the GUI
%         SetGlobalPairCopulaParameterBounds. Furthermore, the function
%         PairCopulaParameterBounds can be used to obtain the current
%         bounds.
%
%
% Usage:
%               CheckParameter(family,theta)
%
%
% Inputs:
%       families  = The copula family.
%       theta     = The parameter vector for the pair-copula.
%
%
% Author: Malte Kurz

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
 
k = 1;
for i = 1:NumbPCs
    switch families(i)
        case {3,4,5,6,12,16,17,19}
            [lb,ub] = PairCopulaParameterBounds(families(i));
            for j=1:2
                if sum(theta(k) < lb(j)) > 0 || sum(theta(k) > ub(j)) > 0
                    throwAsCaller(MException('',['The ' num2str(j) '. parameter of copula family number ' num2str(families(i)) ' has to lie between ' num2str(lb(j)) ' and ' num2str(ub(j)) '.']))
                end
                k=k+1;
            end
            
        case {18}
            [lb,ub] = PairCopulaParameterBounds(families(i));
            for j=1:3
                if sum(theta(k) < lb(j)) > 0 || sum(theta(k) > ub(j)) > 0
                    throwAsCaller(MException('',['The ' num2str(j) '. parameter of copula family number ' num2str(families(i)) ' has to lie between ' num2str(lb(j)) ' and ' num2str(ub(j)) '.']))
                end
                k=k+1;
            end
            
        case {0}
            
        otherwise
            [lb,ub] = PairCopulaParameterBounds(families(i));
            if sum(theta(k) < lb) > 0 || sum(theta(k) > ub) > 0
                throwAsCaller(MException('',['The 1. parameter of copula family number ' num2str(families(i)) ' has to lie between ' num2str(lb) ' and ' num2str(ub) '.']))
            end
            k=k+1;
            
    end
end

end
