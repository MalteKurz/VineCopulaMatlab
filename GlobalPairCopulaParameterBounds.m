function [varargout] = GlobalPairCopulaParameterBounds(varargin)
%GLOBALPAIRCOPULAPARAMETERBOUNDS Setting the global (i.e., for the whole VineCopulaMatlab toolbox) pair-copula parameter bounds, which are used for estimation and consistency checks
% Purpose
%        The function can be used to set the global (i.e., for the whole
%        VineCopulaMatlab toolbox) pair-copula parameter bounds, which are
%        used for estimation and consistency checks. Possible copula
%        families (the default values for the parameter bounds are given in
%        parenthesis, where the ordering is as follows (lb1, ub1, lb2, ub2, lb3, ub3)):
%           0   Indep
%           1   AMH (-1, 1)
%           2   AsymFGM (0, 1)
%           3   BB1 (0, 6, 1, 6)
%           4   BB6 (1, 6, 1, 6)
%           5   BB7 (1, 6, 0.001, 6)
%           6   BB8 (1, 6, 0.001, 1)
%           7   Clayton (0, Inf)
%           8   FGM (-1, 1)
%           9   Frank (-30, 30)
%           10  Gaussian (-0.999, 0.999)
%           11  Gumbel (1, Inf)
%           12  IteratedFGM (-1, 1, -1, 1)
%           13  Joe (1, Inf)
%           14  PartialFrank (0, 30)
%           15  Plackett (0.001, Inf)
%           16  Tawn1 (1.001, 20, 0.001, 0.999)
%           17  Tawn2 (1.001, 20, 0.001, 0.999)
%           18  Tawn (1.0001, 20, 0.001, 0.999, 0.001, 0.999)
%           19  t (-0.999, 0.999, 1, 30)
%
%
% Usage
%               GlobalPairCopulaParameterBounds
%
%
%
% Author: Malte Kurz
mlock

persistent bounds;

bounds_extreme = nan(20,6);
bounds_extreme(2,1:2) = [-1 1]; %AMH
bounds_extreme(3,1:2) = [0 1]; %AsymFGM
bounds_extreme(4,1:4) = [0 Inf 1 Inf]; %BB1
bounds_extreme(5,1:4) = [1 Inf 1 Inf]; %BB6
bounds_extreme(6,1:4) = [1 Inf 0 Inf]; %BB7
bounds_extreme(7,1:4) = [1 Inf 0 1]; %BB8
bounds_extreme(8,1:2) = [0 Inf]; %Clayton
bounds_extreme(9,1:2) = [-1 1]; %FGM
bounds_extreme(10,1:2) = [-Inf Inf]; %Frank
bounds_extreme(11,1:2) = [-1 1]; %Gaussian
bounds_extreme(12,1:2) = [1 Inf]; %Gumbel
bounds_extreme(13,1:4) = [-1 1 -1 1]; %IteratedFGM
bounds_extreme(14,1:2) = [1 Inf]; %Joe
bounds_extreme(15,1:2) = [0 Inf]; %PartialFrank
bounds_extreme(16,1:2) = [0 Inf]; %Plackett
bounds_extreme(17,1:4) = [1 Inf 0 1]; %Tawn1
bounds_extreme(18,1:4) = [1 Inf 0 1]; %Tawn2
bounds_extreme(19,1:6) = [1 Inf 0 1 0 1]; %Tawn
bounds_extreme(20,1:4) = [-1 1 1 Inf]; %t

bounds_extreme = array2table(bounds_extreme,...
    'VariableNames',{'lb1','ub1','lb2','ub2','lb3','ub3'},...
    'RowNames',{'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'});

bounds_default = nan(20,6);
bounds_default(2,1:2) = [-1 1]; %AMH
bounds_default(3,1:2) = [0 1]; %AsymFGM
bounds_default(4,1:4) = [0 6 1 6]; %BB1
bounds_default(5,1:4) = [1 6 1 6]; %BB6
bounds_default(6,1:4) = [1 6 0.001 6]; %BB7
bounds_default(7,1:4) = [1 6 0.001 1]; %BB8
bounds_default(8,1:2) = [0 Inf]; %Clayton
bounds_default(9,1:2) = [-1 1]; %FGM
bounds_default(10,1:2) = [-30 30]; %Frank
bounds_default(11,1:2) = [-0.999 0.999]; %Gaussian
bounds_default(12,1:2) = [1 Inf]; %Gumbel
bounds_default(13,1:4) = [-1 1 -1 1]; %IteratedFGM
bounds_default(14,1:2) = [1 Inf]; %Joe
bounds_default(15,1:2) = [0 30]; %PartialFrank
bounds_default(16,1:2) = [0.001 Inf]; %Plackett
bounds_default(17,1:4) = [1.001 20 0.001 0.999]; %Tawn1
bounds_default(18,1:4) = [1.001 20 0.001 0.999]; %Tawn2
bounds_default(19,1:6) = [1.001 20 0.001 0.999 0.001 0.999]; %Tawn
bounds_default(20,1:4) = [-0.999 0.999 1 30]; %t

bounds_default = array2table(bounds_default,...
    'VariableNames',{'lb1','ub1','lb2','ub2','lb3','ub3'},...
    'RowNames',{'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'});

if isempty(bounds)
    bounds = bounds_default;
end

if nargin==1 % The case where the bounds should be updated
    % If a table is provided, the families are checked first.
    if istable(varargin{1})
        UnknownFamilies = setdiff(varargin{1}.Properties.RowNames,bounds.Properties.RowNames);
        if not(isempty(UnknownFamilies))
            error(['Unknown copula familie(s) ' strjoin(UnknownFamilies)])
        end
        bounds(varargin{1}.Properties.RowNames,:) = varargin{1};
        
        % Check whether the parameter bounds are set correct
        if any(any(table2array(bounds(:,[1,3,5]))<table2array(bounds_extreme(:,[1,3,5])))) || any(any(table2array(bounds(:,[2,4,6]))>table2array(bounds_extreme(:,[2,4,6]))))
            error(['Invalid parameter bound(s) for the copula familie(s) ' strjoin(bounds.Properties.RowNames(any([any(table2array(bounds(:,[1,3,5]))<table2array(bounds_extreme(:,[1,3,5])),2),any(table2array(bounds(:,[2,4,6]))>table2array(bounds_extreme(:,[2,4,6])),2)],2)))])
        end
    else
        % FIXME Place code for matrix input here
    end
    
    varargout{1} = table2array(bounds);
else
    varargout{1} = table2array(bounds);
end



end
