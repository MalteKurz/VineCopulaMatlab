function ParamHat = PairCopulaFit(family,u1,u2,varargin)
%PAIRCOPULAFIT Estimating pair-copulas
% Purpose
%        The function returns the ML-estimates for two-dimensional (i.e.,
%        pair-) copulas. Possible copula families:
%           0   Indep
%           1   AMH
%           2   AsymFGM
%           3   BB1
%           4   BB6
%           5   BB7
%           6   BB8
%           7   Clayton
%           8   FGM
%           9   Frank
%           10  Gaussian
%           11  Gumbel
%           12  IteratedFGM
%           13  Joe
%           14  PartialFrank
%           15  Plackett
%           16  Tawn1
%           17  Tawn2
%           18  Tawn
%           19  t
%
%
% Usage
%               ParamHat = PairCopulaFit(family,u1,u2)
%           Rotated pair-copulas
%               ParamHat = PairCopulaFit(family,u1,u2,rotation)
% 
% 
% Inputs
%       family    = The copula family.
%       u1        = A (n x 1) dimensional vector of values lying in [0,1].
%       u2        = A (n x 1) dimensional vector of values lying in [0,1].
%       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
%                   rotation is achieved by letting the rotation argument
%                   empty or by choosing 0 rotation.
%
%
% Outputs
%       ParamHat  = The ML-estimate for the parameter vector.
%
%
%
% Author: Malte Kurz

% Check the (Copula-)data input.
CheckData([u1,u2])

if not(isnumeric(family))
    families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
    if sum(strcmpi(family,families))
        family = find(strcmp(family,families))-1;
    else
        error(['The copula family ' family ' is not implemented'])
    end
end
        

% Checking for rotated copulas
if nargin == 4 && not(isempty(varargin{1}))
    if not(isnumeric(varargin{1}))
        Rotations = {'r90','r180','r270'};
        varargin{1} = find(strcmp(varargin{1},Rotations)).*90;
    end
    ParamHat = PCFit(family,u1,u2,varargin{1});
else
    ParamHat = PCFit(family,u1,u2);
end

end
