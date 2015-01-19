function U = PairCopulaHfun(family,u1,u2,varargin)
%PAIRCOPULAHFUN Computing the h-function of a pair-copula
% Purpose
%        The function computes the h-function, i.e., it evaluates the
%        conditional distribution function of U1|U2 at u1, given U2=u2.
%        Possible copula families:
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
%               U = PairCopulaHfun(family,u1,u2,theta)
%           Rotated pair-copulas
%               U = PairCopulaHfun(family,u1,u2,theta,rotation)
%
%
% Inputs
%       family    = The copula family.
%       u1        = A (n x 1) dimensional vector of values lying in [0,1].
%       u2        = A (n x 1) dimensional vector of values lying in [0,1].
%       theta     = The parameter vector for the pair-copula.
%       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
%                   rotation is achieved by letting the rotation argument
%                   empty or by choosing 0 rotation.
%
%
% Outputs
%      U          = n dimensional vector of probabilities (sometimes 
%                   called pseudo-observations), which are lying strictly
%                   in [0,1].
%
%
%
% Author: Malte Kurz

% Checking the dimensions
if not(size(u1,2) == 1)
    error('The u1 vector has to be of dimension n x 1.')
end
if not(size(u1,1) == size(u2,1))
    error('The number of observations (rows) of the vectors u1 and u2 have to be equal.')
end

% Check the (Copula-)data input.
CheckData([u1,u2])

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

if family == 0
    U=PCHfun(0,u1,u2,0);
elseif nargin == 5 && not(isempty(varargin{2}))
    if not(isnumeric(varargin{2}))
        Rotations = {'r90','r180','r270'};
        varargin{2} = find(strcmp(varargin{2},Rotations)).*90;
    end
    U=PCHfun(family,u1,u2,varargin{1},varargin{2});
else
    U=PCHfun(family,u1,u2,varargin{1});
end

end
