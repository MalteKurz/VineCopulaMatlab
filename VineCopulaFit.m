function [ParamHat,varargout] = VineCopulaFit(type,families,d,u,varargin)
%VINECOPULAFIT Estimating vine copulas
% Purpose
%        The function computes ML-estimates for the parameters of a C-Vine
%        or D-Vine copula. Therefore, first starting values for the joint
%        estimation are obtained by iteratively estimating the pair-copulas
%        in the first trees and using those estimates to obtain the 
%        arguments for the copulas in the second tree. Then the pair-
%        copulas in the second tree are estimated and so on. These
%        sequentially estimated parameters from the sequential procedure
%        are then used to obtain the ML-estimates, by minimizing the
%        overall negative  log-likelihood of the whole C-Vine (or D-Vine)
%        numerically. Possible Vine copula types:
%           0   C-Vine
%           1   D-Vine
%
%
% Usage
%           Simplified standard C-Vine or D-Vine copula
%               ParamHat = VineCopulaFit('C-Vine',families,d,u)
%               [ParamHat, MaxLogLikes] = VineCopulaFit('C-Vine',families,d,u)
%               [ParamHat, MaxLogLikes, theta0] = VineCopulaFit('C-Vine',families,d,u)
%           Simplified standard C-Vine or D-Vine copula with rotated copulas
%               ParamHat = VineCopulaFit('C-Vine',families,d,u,rotation)
%               [ParamHat, MaxLogLikes] = VineCopulaFit('C-Vine',families,d,u,rotation)
%               [ParamHat, MaxLogLikes, theta0] = VineCopulaFit('C-Vine',families,d,u,rotation)
%           Simplified standard C-Vine or D-Vine copula (specified estimation method
%           (i.e., joint or sequential estimation)
%               ParamHat = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod)
%               [ParamHat, MaxLogLikes] = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod)
%               [ParamHat, MaxLogLikes, theta0] = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod)
%           Truncated simplified standard C-Vine or D-Vine copula
%               ParamHat = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod,CutOffTree)
%               [ParamHat, MaxLogLikes] = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod,CutOffTree)
%               [ParamHat, MaxLogLikes, theta0] = VineCopulaFit('C-Vine',families,d,u,rotation,EstMethod,CutOffTree)
%
%
% Inputs
%       type            = The vine copula type.
%       families        = A vector of the pair-copula families, which
%                         are part of the PCC. The vector has to have
%                         the length (d-1)*d/2. The first d-1 entries are
%                         the copula families in the first tree and the
%                         next d-2 entries are the copula families in the
%                         second tree and so on. That means, for d=4 the
%                         array should look similar to this {'Frank',
%                         'Frank', 'Frank', 'AMH', 'AMH', 'Clayton'}, which
%                         is the special case where all copulas in the
%                         first tree are Frank copulas, all copulas in the
%                         second tree are AMH copulas and all copulas in
%                         the third tree are Clayton copulas.
%                         The order of the families is:
%                         * Exemplarily for the four-dimensional C-Vine):
%                           C12, C13, C14, C23|1, C24|1, C34|12
%                         * Exemplarily for the four-dimensional D-Vine):
%                           C12, C23, C34, C13|2, C24|3, C14|23
%                         Note: If families is a simple string/character,
%                         e.g., 'Clayton', then all pair-copulas are
%                         specified to be from this copula family.
%       d               = The dimension of the C- or D-Vine.
%       u               = A (n x d) dimensional vector of values lying in
%                         [0,1] (the observations).
%       rotation        = A vector of the same dimension as families in
%                         which one can specify rotation levels.
%       EstMethod       = The estimation method must be either 'joint' or
%                         'sequential'. If it is not explicitly given, a
%                         joint estimation is performed (default).
%       CutOffTree      = The CutOffTree (or also called truncation level)
%                         can be used to set all pair-copulas from the
%                         (CutOffTree + 1)-th tree on to independence
%                         copulas (i.e., ignore them in the joint
%                         estimation).
%
%
% Outputs
%       ParamHat        = The ML-estimates of the parameters for the
%                         (d-1)*d/2 pair-copulas. These estimates are given
%                         in the same order as the families vector, but
%                         in a row-vector. If a pair-copula is an
%                         independence copula, then there is no estimate
%                         given. Furthermore, if a pair-copula has two or 
%                         more parameters, the estimates are given in same 
%                         order as they have to be provided if the pair-
%                         copula is considered only. For example, for a t-
%                         copula, the first estimate is for the parameter
%                         rho and the second one for the degrees of freedom
%                         parameter nu.
%       MaxLogLikes     = The first entry is the value of copula-log-
%                         likelihood evaluated at the sequentially estimated
%                         ML-estimates. The second value is the value of
%                         the copula-log-likelihood evaluated for the
%                         joint ML-estimates.
%       theta0          = The vector of ML-estimates, which are obtained by
%                         using the sequential estimation approach. Note
%                         that these estimates are also used as starting
%                         point for the global maximum likelihood
%                         estimation.
%
%
%
% Author: Malte Kurz

% Check the (Copula-)data input.
CheckData(u)

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

if nargin == 7
    CutOffTree = varargin{3};
    K = (d-1) -min(d-1,CutOffTree);
    NumbPCs = d*(d-1)/2 - K*(K+1)/2;
    if not(isempty(varargin{2}))
        EstMethod = varargin{2};
    else
        EstMethod = 'joint';
    end
    families(1,NumbPCs+1:d*(d-1)/2) = 0; % Set the copulas in the upper trees to independence copulas.
    if not(isempty(varargin{1}))
        if not(isnumeric(varargin{1}))
            Rotations = {'r90','r180','r270'};
            for i=1:length(varargin{1})
                if isempty(varargin{1}{i})
                    varargin{1}{i} = 0;
                else
                    varargin{1}{i} = find(strcmp(varargin{1}(i),Rotations)).*90;
                end
            end
            varargin{1} = cell2mat(varargin{1});
        end
        varargin{1}(NumbPCs+1:d*(d-1)/2) = 0;
    else
        varargin{1} = zeros(1,d*(d-1)/2);
    end
elseif nargin == 6
    if not(isempty(varargin{2}))
        EstMethod = varargin{2};
    else
        EstMethod = 'joint';
    end
    NumbPCs = d*(d-1)/2;
    if isempty(varargin{1})
        varargin{1} = zeros(1,d*(d-1)/2);
    end
else
    EstMethod = 'joint';
    NumbPCs = d*(d-1)/2;
end

if ~(strcmp(EstMethod,'joint') || strcmp(EstMethod,'sequential'))
    error('The variable EstMethod must be either ''joint'' or ''sequential''.')
end

% Checking for the right number of copula families
if length(families)==1
    families = repmat(families,1,d*(d-1)/2);
elseif (length(families) ~= d*(d-1)/2) && (length(families) ~= NumbPCs);
    error('Wrong number of copula families')
end

if not(isnumeric(varargin{1}))
    Rotations = {'r90','r180','r270'};
    for i=1:length(varargin{1})
        if isempty(varargin{1}{i})
            varargin{1}{i} = 0;
        else
            varargin{1}{i} = find(strcmp(varargin{1}(i),Rotations)).*90;
        end
    end
    varargin{1} = cell2mat(varargin{1});
end

if strcmp(EstMethod,'joint')
    [ParamHat,MaxLogLikes,theta0] = VineFit(type,u,families,varargin{1});
else
    [theta0,MaxLogLikes] = VineFitSeq(type,u,families,varargin{1});
    ParamHat = theta0;
end

varargout{1} = MaxLogLikes;
varargout{2} = theta0;
if nargin == 7
    varargout{3} = [families;varargin{1}];
end

end
