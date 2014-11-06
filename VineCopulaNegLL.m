function CLL = VineCopulaNegLL(type,u,families,thetas,d,varargin)
%VINECOPULANEGLL Computing the negative log-likelihood of a vine copula
% Purpose
%        The function computes the value of the log-likelihood for a given 
%        matrix of observations u, which have to lie in the d-dimensional
%        unit cube. Possible Vine copula types:
%           0   C-Vine
%           1   D-Vine
%
%
% Usage
%           Simplified standard C-Vine or D-Vine copula
%               CLL = CopulaNegLL('C-Vine',u,families,thetas,d)
%               CLL = CopulaNegLL('D-Vine',u,families,thetas,d)
%           Simplified standard C-Vine or D-Vine copula with rotated copulas
%               CLL = CopulaNegLL('C-Vine',u,families,thetas,d,rotation)
%               CLL = CopulaNegLL('D-Vine',u,families,thetas,d,rotation)
%           Truncated simplified standard C-Vine or D-Vine copula
%               CLL = CopulaNegLL('C-Vine',u,families,thetas,d,rotation,CutOffTree)
%               CLL = CopulaNegLL('D-Vine',u,families,thetas,d,rotation,CutOffTree)
%
%
% Inputs
%       type        = The vine copula type.
%       u           = A (n x d) dimensional vector of values lying in [0,1]
%                     (the observations).
%       families    = A vector of the pair-copula families, which
%                     are part of the PCC. The vector has to have
%                     the length (d-1)*d/2. The first d-1 entries are
%                     the copula families in the first tree and the
%                     next d-2 entries are the copula families in the
%                     second tree and so on. That means, for d=4 the
%                     array should look similar to this {'Frank',
%                     'Frank', 'Frank', 'AMH', 'AMH', 'Clayton'}, which
%                     is the special case where all copulas in the
%                     first tree are Frank copulas, all copulas in the
%                     second tree are AMH copulas and all copulas in
%                     the third tree are Clayton copulas.
%                     The order of the families is:
%                     * Exemplarily for the four-dimensional C-Vine):
%                           C12, C13, C14, C23|1, C24|1, C34|12
%                     * Exemplarily for the four-dimensional D-Vine):
%                           C12, C23, C34, C13|2, C24|3, C14|23
%                     Note: If families is a simple string/character,
%                     e.g., 'Clayton', then all pair-copulas are
%                     specified to be from this copula family.
%       thetas      = The values of the parameters for the (d-1)*d/2 pair-
%                     copulas. These parameter values have to be given in
%                     the same order as the families vector, but in a
%                     row-vector. If a pair-copula is a independence
%                     copula, then there is no parameter needed.
%                     Furthermore, if a pair-copula has two or more
%                     parameters, the parameters have to be given in same 
%                     order as they have to be provided if the pair-copula
%                     is considered only. For example for a t-copula, the
%                     first parameter is rho and the second parameter is
%                     the degrees of freedom parameter nu.
%                     Note: If thetas is a skalar then all parameters
%                     (i.e., the parameters of all pair-copulas) are set to
%                     the same value (i.e. to the one given).
%       d           = The dimension of the C- or D-Vine.
%       rotation    = A vector of the same dimension as families in
%                     which one can specify rotation levels.
%       CutOffTree  = The CutOffTree (or also called truncation level) can
%                     be used to set all pair-copulas from the (CutOffTree
%                     + 1)-th tree on to independence copulas (i.e., ignore
%                     them in the joint estimation).
%
% Outputs
%      CLL          = The value of the negative log-likelihood for
%                     the data matrix u.
%
%
%
% Author: Malte Kurz

% Check the (Copula-)data input.
CheckData(u)
 
% Check the copula parameters.
CheckParameters(families,thetas)

narginchk(5, 7)

if not(isnumeric(type))
    type = ceil(find(strcmp(type,{'CVine','C-Vine','DVine','D-Vine'}))/2)-1;
end

if not(type == 0 || type == 1)
    error('Invalid vine copula type')
end

if not(isnumeric(families))
    Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
    for i=1:length(families)
        if sum(strcmpi(families(i),Families))
            families{i} = find(strcmp(families(i),Families))-1;
        else
            error(['The copula family ' families{i} ' is not implemented'])
        end
    end
    families = cell2mat(families);
end
    
if nargin == 6
    NumbPCs = d*(d-1)/2;
    if not(isempty(varargin{1}))
        rotation = varargin{1};
    else
        rotation = zeros(1,NumbPCs);
    end
    CutOffTree = d-1;
elseif nargin == 7
    if not(isempty(varargin{1}))
        rotation = varargin{1};
    else
        rotation = zeros(1,d*(d-1)/2);
    end
    if not(isempty(varargin{2}))
        CutOffTree = varargin{2}; % The last tree, which should be included into the likelihood.
    else
        CutOffTree = d-1;
    end
    J = min(d-1,CutOffTree);
    K = (d-1) - J;
    NumbPCs = d*(d-1)/2 - K*(K+1)/2;
else
    NumbPCs = d*(d-1)/2;
    rotation = zeros(1,NumbPCs);
    CutOffTree = d-1;
end

if not(isnumeric(rotation))
    Rotations = {'r90','r180','r270'};
    for i=1:length(rotation)
        if isempty(rotation{i})
            rotation{i} = 0;
        else
            rotation{i} = find(strcmp(rotation(i),Rotations)).*90;
        end
    end
    rotation = cell2mat(rotation);
end

% Checking for the right number of copula families
if length(families)==1
    families = repmat(families,1,d*(d-1)/2);
elseif (length(families) ~= d*(d-1)/2) && (length(families) ~= NumbPCs);
    error('Wrong number of copula families')
end

if length(thetas) == 1
    thetas = repmat(thetas,1,NumbPCs + sum(families(1:NumbPCs)==3) + sum(families(1:NumbPCs)==4) + sum(families(1:NumbPCs)==5) + sum(families(1:NumbPCs)==6) + sum(families(1:NumbPCs)==12) + sum(families(1:NumbPCs)==16) + sum(families(1:NumbPCs)==17) + sum(families(1:NumbPCs)==19)+ 2*sum(families(1:NumbPCs)==18) - sum(families(1:NumbPCs)==0));
elseif length(thetas) ~= NumbPCs+ sum(families(1:NumbPCs)==3) + sum(families(1:NumbPCs)==4) + sum(families(1:NumbPCs)==5) + sum(families(1:NumbPCs)==6) + sum(families(1:NumbPCs)==12) + sum(families(1:NumbPCs)==16) + sum(families(1:NumbPCs)==17) + sum(families(1:NumbPCs)==19)+ 2*sum(families(1:NumbPCs)==18) - sum(families(1:NumbPCs)==0)
    error('Wrong number of copula parameters')
end

CLL = VineNegLL(type,u,families,thetas,rotation,CutOffTree);

end
