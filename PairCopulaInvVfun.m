function U = PairCopulaInvVfun(family,u1,u2,varargin)
%PAIRCOPULAINVVFUN Computing the inverse of the v-function of a pair-copula
% Purpose
%        The function computes the inverse of the v-function, (i.e., the 
%        inverse of the conditional distribution function of U2|U1 in the
%        first argument, given U1=u1, evaluated at u2. Possible copula
%        families:
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
%               U = PairCopulaInvVfun(family,u1,u2,theta)
%           Rotated pair-copulas
%               U = PairCopulaInvVfun(family,u1,u2,theta,rotation)
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
%      p          = The value of the inverse v-function evaluated at u2,
%                   given U1 = u1.
%
%
%
% Author: Malte Kurz

% Checking the dimensions
if not(size(u2,2) == 1)
    error('The u2 vector has to be of dimension n x 1.')
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
    U=PCInvVfun(0,u1,u2,0);
elseif nargin == 5 && not(isempty(varargin{2}))
    if not(isnumeric(varargin{2}))
        Rotations = {'r90','r180','r270'};
        varargin{2} = find(strcmp(varargin{2},Rotations)).*90;
    end
    U=PCInvVfun(family,u1,u2,varargin{1},varargin{2});
else
    U=PCInvVfun(family,u1,u2,varargin{1});
end

end
