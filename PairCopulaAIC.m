function [AIC,ParamHat] = PairCopulaAIC(family,u1,u2,varargin)
%PAIRCOPULAAIC Computing the AIC of a pair-copula
% Purpose
%        The function computes the value of the AIC for a copula for a
%        given matrix of observations u, which have to lie in the
%        2-dimensional unit cube, evaluated at the ML estimates. Possible
%        pair-copula families:
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
%        [AIC,ParamHat] = PairCopulaAIC(family,u1,u2)
%      Rotated pair-copulas
%        [AIC,ParamHat] = PairCopulaAIC(family,u1,u2,roatation)
%
%
% Inputs
%       family          = The copula family.
%       u1              = A (n x 1) dimensional vector of values lying in [0,1].
%       u2              = A (n x 1) dimensional vector of values lying in [0,1].
%       rotation        = The degree of rotation, i.e., either 90, 180 or
%                         270. No rotation is achieved by letting the
%                         rotation argument empty or by choosing 0
%                         rotation.
%
%
% Outputs
%      AIC              = The value of the AIC evaluated at the ML-estimator.
%      ParamHat         = The ML-estimate for the parameter vector.
%
%
%
% Author: Malte Kurz

% Checking the dimensions
if not(size(u1,2) == 1 && size(u2,2) == 1)
    error('The u1 and u2 vector have to be of dimension n x 1.')
end
if not(size(u1,1) == size(u2,1))
    error('The number of observations (rows) of the vectors u1 and u2 have to be equal.')
end

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

% Obtain parameter bounds
bounds = GlobalPairCopulaParameterBounds';

if nargin == 4 && not(isempty(varargin{1}))
    if not(isnumeric(varargin{1}))
        Rotations = {'r90','r180','r270'};
        varargin{1} = find(strcmp(varargin{1},Rotations)).*90;
    end
    [AIC,ParamHat] = VineCopulaMatlab(1,bounds,family,u1,u2,varargin{1});
else
    [AIC,ParamHat] = VineCopulaMatlab(1,bounds,family,u1,u2);
end
        
end
