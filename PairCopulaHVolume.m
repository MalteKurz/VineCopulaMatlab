function P = PairCopulaHVolume(family,a,b,varargin)
%PAIRCOPULAHVOLUME Computing the h-volume of a copula
% Purpose
%        The function computes the probability for a random vector, being
%        distributed according to a specific copula, to lie in a
%        hyperrectangle. The hyperrectangle is defined by the cartesean
%        product of the intervals specified by the lower bounds a and upper
%        bounds b. Possible pair-copula families:
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
%               P = CopulaHVolume(family,a,b,theta)
%           Rotated pair-copulas
%               P = CopulaHVolume(family,a,b,theta,rotation)
%
%
% Inputs
%       family    = The copula family.
%       a         = A 2-dimensional vector of lower bounds for 2 intervals,
%                   defining a 2-dimensional hyperrectangle.
%       b         = A 2-dimensional vector of upper bounds for 2 intervals,
%                   defining a 2-dimensional hyperrectangle.
%       theta     = The parameter vector for the pair-copula.
%       rotation  = The degree of rotation, i.e., either 90, 180 or 270. No
%                   rotation is achieved by letting the rotation argument
%                   empty or by choosing 0 rotation.
%
%
% Outputs
%      P          = The probability for a random variable, which is
%                   distributed according to the specified pair-copula, to
%                   lie in the hyperrectangle defined by the lower and
%                   upper bounds vectors a and b.
%
%
%
% Author: Malte Kurz

% Check for correct inputs
if sum(a>=b) || sum(b>1) || sum(a<0)
    error('Wrong input for the bounds.')
end

% Check the (Copula-)data input.
CheckData([a,b])

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
    C = @(U) VineCopulaMatlab(2,0,U(:,1),U(:,2));
elseif nargin == 5 && not(isempty(varargin{2}))
    if not(isnumeric(varargin{2}))
        Rotations = {'r90','r180','r270'};
        varargin{2} = find(strcmp(varargin{2},Rotations)).*90;
    end
    C = @(U) VineCopulaMatlab(2,family,U(:,1),U(:,2),varargin{1},varargin{2});
else
    C = @(U) VineCopulaMatlab(2,family,U(:,1),U(:,2),varargin{1});
end

dim = length(a);
% Transfer the vectors of bounds into row vectors.
a = a(:)';
b = b(:)';

if dim == 2
    upos = [b ; a];
    uneg = [b(1) a(2); a(1) b(2)];
%elseif dim==3
%    upos = [b ; b(1) a(2:3); a(1) b(2) a(3); a(1:2) b(3)];
%    uneg = [b(1:2) a(3);b(1) a(2) b(3); a(1) b(2:3); a];
else
    error('Only the two-dimensional case is implemented.')
end

P = sum(C(upos))-sum(C(uneg));

end
