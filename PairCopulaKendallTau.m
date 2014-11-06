function tau = PairCopulaKendallTau(family,theta)
%PAIRCOPULAKENDALLTAU Transforming the parameter of a pair-copula into the value of Kendall's tau
% Purpose
%        The function computes the value of Kendall's tau given the
%        parameter of the specific two-dimensional copula. Possible copula
%        families:
%           1   AMH
%           7   Clayton
%           8   FGM
%           9   Frank
%           10  Gaussian
%           11  Gumbel
%           19  t
%
%
% Usage
%               tau = CopulaKendallTau(family,theta)
%
%
% Inputs
%       family      = The copula family.
%       theta       = The parameter vector for the pair-copula.
%
%
% Outputs
%      tau          = The value of Kendall's tau.
%
%
%
% Author: Malte Kurz

if not(isnumeric(family))
    families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
    if sum(strcmpi(family,families))
        family = find(strcmp(family,families))-1;
    else
        error(['The copula family ' family ' is not implemented'])
    end
end

% Check the copula parameters.
if family==19 && length(theta) == 2
    % Set the degrees of freedom parameter to an allowed value to pass the
    % CheckParameters test.
    theta(2) = 3;
end
CheckParameters(family,theta)

switch family
    case 1 % AMH
        tau = (3.*theta-2)./(3.*theta)-(2.*(1-theta).^2.*log(1-theta))./(3.*theta.^2);
        
    case 7 % Clayton
        tau = theta./(theta+2);
        
    case 8 % FGM
        tau = 2/9*theta;
        
    case 9 % Frank
        DebyeIntegrand = @(t) t./expm1(t);
        DebyeFun = zeros(length(theta),1);
        for i = 1:length(theta)
            if theta(i) == 0
                % Special case for the Frank copula, which is the
                % independence copula if theta=0
                DebyeFun(i,1) = 0;
                theta(i) = 4;
            else
                DebyeFun(i,1) = 1./theta(i)*quad(DebyeIntegrand,0,theta(i));
            end
        end
        tau = 1-4./theta.*(1-DebyeFun);
    
    case 10 % Gaussian
        tau = 6/pi.*asin(0.5.*theta);
        
    case 11 % Gumbel
        tau = 1-1./theta;
        
    case 19 % t
        tau = 2/pi.*asin(theta);
        
    otherwise
        error(['The copula family ' families{family+1} ' is not implemented.'])
        
end

end
