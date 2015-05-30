function U = PairCopulaRand(family,N,varargin)
%PAIRCOPULARAND Generating (pseudo-)random variables from a pair-copula
% Purpose
%        The function draws N pseudo-random 2-dimensional tuples from a
%        pair-copula. Possible copula families:
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
%               U = CopulaRand(family,N,theta)
%           Rotated pair-copulas
%               U = CopulaRand(family,N,theta,rotation)
%
%
% Inputs
%       family    = The copula family.
%       N         = The number of pseudo-random tuples to be drawn.
%       theta     = The parameter vector for the pair-copula.
%       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
%                   rotation is achieved by letting the rotation argument
%                   empty or by choosing 0 rotation.
%
%
% Outputs
%      U          = (N x 2) Matrix of simulated tuples from the specified
%                   pair-copula.
%
%
%
% Author: Malte Kurz

% Check the copula parameters.
CheckParameters(family,varargin{1})

if not(isnumeric(family))
    families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
    if sum(strcmpi(family,families))
        family = find(strcmp(family,families))-1;
    else
        error(['The copula family ' family ' is not implemented'])
    end
end

U = zeros(N,2);

if family == 0
    [U(:,1),U(:,2)] = VineCopulaMatlab(10,0,N);
elseif nargin == 4 && not(isempty(varargin{2}))
    if not(isnumeric(varargin{2}))
        Rotations = {'r90','r180','r270'};
        varargin{2} = find(strcmp(varargin{2},Rotations)).*90;
    end
    [U(:,1),U(:,2)] = VineCopulaMatlab(10,family,N,varargin{1},varargin{2});
else
    [U(:,1),U(:,2)]= VineCopulaMatlab(10,family,N,varargin{1});
end

end
