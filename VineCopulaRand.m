function U = VineCopulaRand(type,N,d,families,thetas,varargin)
%VINECOPULARAND Generating (pseudo-)random variables from a copula or vine copula
% Purpose
%        The function draws N pseudo-random d-dimensional tuples from a
%        C-Vine or D-Vine copula. The vine copulas can be of the simplified
%        and non-simplified form. In the case of a non-simplified Vine
%        copula, the parameters of the conditional copulas have to be given
%        as functionals of the variables on which the conditioning is done.
%
%
% Usage
%           Standard C-Vine copula (simplified & non-simplified, where the
%           parameters of the conditional copulas have to be given as
%           functionals of the variables on which the conditioning is done)
%             simplified
%               U = CopulaRand('C-Vine',N,d,families,thetas,ones(1,(d-1)*(d-2)/2))
%             non-simplified
%               U = CopulaRand('C-Vine',N,d,families,thetas,simplified,condparameterfunctionals)
%           Standard D-Vine copula (simplified & non-simplified, where the
%           parameters of the conditional copulas have to be given as
%           functionals of the variables on which the conditioning is done)
%             simplified
%               U = CopulaRand('D-Vine',N,d,families,thetas,ones(1,(d-1)*(d-2)/2))
%             non-simplified
%               U = CopulaRand('D-Vine',N,d,families,thetas,simplified,condparameterfunctionals)
%
%
% Inputs
%       type      = The vine copula type.
%       N         = The number of pseudo-random tuples to be drawn.
%       d         = The dimension of the C- or D-Vine.
%       rotation    = A vector of the same dimension as families in
%                     which one can specify rotation levels.
%       families  = A cell array of the pair-copula families, which are
%                   part of the PCC. The cell array has to have the length
%                   (d-1)*d/2. The first d-1 entries are the copula
%                   families in the first tree and the next d-2 entries are
%                   the copula families in the second tree and so on. That
%                   means, for d=4 the array should look similar to this
%                   {'Frank', 'Frank', 'Frank', 'AMH', 'AMH', 'Clayton'},
%                   which is the special case where all copulas in the
%                   second tree are AMH copulas and all copulas in the
%                   third tree are Clayton copulas. The order of the
%                   families is (exemplarily for the case d=4): C12, C13,
%                   C14, C23|1, C24|1, C34|12.
%                   Note: If families is a simple string/character, e.g.,
%                   'Clayton', then all pair-copulas are specified to be
%                   from this copula family.
%       thetas    = The values of the parameters for the (d-1)*d/2 pair-
%                   copulas. These parameter values have to be given in the
%                   same order as the families vector, but in a row-
%                   vector. If a pair-copula is a independence copula, then
%                   there is no parameter needed. Furthermore, if a pair-
%                   copula has to or more parameters, the parameters have
%                   to be given in same order as they have to be provided
%                   if the pair-copula is considered only. For example for
%                   a t-copula, the first parameter is rho and the second
%                   parameter is the degrees of freedom parameter nu.
%                   Note: If thetas is a skalar then all parameters (i.e.,
%                   the parameters of all pair-copulas) are set to the same
%                   value (i.e. to the one given).
%       simplified
%                 = A vector consisting of (d-1)*(d-2)/2 zeros or ones
%                   (either 1 (for an unconditional bivariate pair-copula)
%                   or 0 for a conditional copula) specifying, whether
%                   the copulas in the second and higher tree are pair-copulas
%                   (unconditional bivariate copulas) or conditional copulas.
%       condparameterfunctionals
%                 = A cell array consisting of (# non-partial (i.e.
%                   conditional) pair-copulas) cell array, which consist of
%                   as many parameter functionals as the copula families
%                   has parameters. Therefore, if one of the parameters is
%                   not varying in the conditioning set, then the function
%                   handle has to give back constant (e.g., @(w) c, where c
%                   is a scalar constant. Functionals for conditional
%                   copulas beeing part of the second tree of the C-Vine
%                   have to be functions of one variable (i.e., they should
%                   be evaluatable for column vectors). Functionals of the
%                   third tree of the C-Vine have to be functionals of two
%                   variables (i.e. they should be evaluatable for matrices
%                   of the size N x 2). And so on, for conditional copulas,
%                   where the conditioning is done on three or more
%                   variables. The ordering of the variables in the
%                   conditioning set corresponds to the ordering of the
%                   root nodes of the different trees (e.g., in the fourth
%                   tree, the first conditional copula is the C_45|123
%                   copula and the conditioning set is three-dimensional,
%                   where the first column of observations corresponds to
%                   observations of the root node in the first tree (i.e.
%                   variable X_1).
%                   Note: Of course it is possible, that the parameter of a
%                   conditional copula only depends on some of the
%                   variables, but then the functions have to be specified
%                   in an appropriate way (e.g., if a conditional copula in
%                   the third tree, where the conditioning is done on two
%                   variables, only depends on the second variable, the
%                   function could look similar to this @(u) 3.*u(:,2),
%                   when the parameter of the conditional copula is exactly
%                   three times the second variable on which the
%                   conditioning is done.
%
%
% Outputs
%      U          = (N x d) Matrix of simulated tuples from the specified
%                   copula, where every row is a simulated d-dimensional
%                   tuple from the d-dimensional C-Vine or D-Vine copula.
%
%
%
% Author: Malte Kurz

% Check the copula parameters.

narginchk(6, 8)

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

simplified = varargin{1};
if sum(simplified==0)==0
CheckParameters(families,thetas)
    if nargin == 7 && not(isempty(varargin{2}))
        error('It is not possible to simulate from a simplified vine-copula model, which contains conditional (non-partial) copulas with functional parameters.')
    end
    if nargin == 8 && not(isempty(varargin{3}))
        rotation = varargin{3};
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
    else
        rotation = zeros(1,d*(d-1)/2);
    end
    % The simplified case is completely simulated in C++
    U=VineCopulaMatlab(105,type,N,d,families,thetas,rotation);
else
    condparameterfunctionals = varargin{2};
    if nargin == 8 && not(isempty(varargin{3}))
        rotation = varargin{3};
    else
        rotation = zeros(1,d*(d-1)/2);
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
    
    
    
    switch type
        case 0
            % In the cell-array InvvfFunctionals for all pair-copulas the
            % inverse v-functions are stored.
            InvvFunctionals = cell(1,d*(d-1)/2);
            
            j = 1;
            k = 1;
            for i = 1:d*(d-1)/2
                if i < d
                    % Here the inverse v-functions for the unconditional
                    % pair-copulas as part of the first tree of the C-Vine are
                    % obtained.
                    switch families(i)
                        case {3,4,5,6,12,16,17,19}
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1),rotation(i));
                            end
                            j = j+2;
                            
                        case {18}
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2),rotation(i));
                            end
                            j = j+3;
                            
                        case {0}
                            InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,0);
                            
                        otherwise
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j),rotation(i));
                            end
                            j = j+1;
                            
                    end
                else
                    % Here the inverse v-functions for the second and higher
                    % trees are obtained, where for each pair-copula it is
                    % distinguished between the pair-copulas, which are partial
                    % copulas and the pair-copulas which are conditional
                    % copulas.
                    if simplified(i-d+1)
                        switch families(i)
                            case {3,4,5,6,12,16,17,19}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1),rotation(i));
                                end
                                j = j+2;
                                
                            case {18}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2),rotation(i));
                                end
                                j = j+3;
                                
                            case {0}
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,0);
                                
                            otherwise
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j),rotation(i));
                                end
                                j = j+1;
                                
                        end
                    else
                        switch families(i)
                            case {3,4,5,6,12,16,17,19}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)]);
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)],rotation(i));
                                end
                                k = k+2;
                                
                            case {18}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)]);
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)],rotation(i));
                                end
                                k = k+3;
                                
                            case {0}
                                InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,0);
                                
                            otherwise
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,condparameterfunctionals{k}(u3));
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,condparameterfunctionals{k}(u3),rotation(i));
                                end
                                k = k+1;
                                
                        end
                    end
                end
            end
            
            % In the indices vector it is stored in which order the inverse
            % h-functions corresponding to the different pair-copulas have to
            % be applied to obtain a pseudo-random sample from the specified
            % C-Vine.
            ind = zeros(1,d*(d-1)/2);
            ind(1:3) = [1 d 2];
            
            if d > 3
                j = 4;
                for i = 3:d-1
                    ind(j:j+i-1) = [ind(j-i+1)+d-i+1,ind(j-i+1:j-1)+1];
                    j = j+i;
                end
            end
            
            % Here the pseudo-random sample from the specified C-Vine is
            % obtained by applying the iterative inversion algorithm for
            % generating pseudo-random variables.
            k = 2;
            W = RandUniform(N,d);
            U(:,1) = W(:,1);
            U(:,2) = InvvFunctionals{1}(W(:,1),W(:,2));
            for i = 3:d
                t = W(:,i);
                for j = i-1:-1:1
                    if ind(k)<d || (simplified(ind(k)-d+1)==1)
                        t = InvvFunctionals{ind(k)}(W(:,j),t);
                    else
                        t = InvvFunctionals{ind(k)}(W(:,j),t,U(:,1:j-1));
                    end
                    k = k+1;
                end
                U(:,i) = t;
            end
            
        case 1
            % In the cell-array InvvfFunctionals for all pair-copulas the
            % inverse v-functions are stored. Furthermore, the cell-array
            % hFunctionals contains all the h-functions;
            InvvFunctionals = cell(1,d*(d-1)/2);
            hFunctionals = cell(1,d*(d-1)/2);
            
            j = 1;
            k = 1;
            for i = 1:d*(d-1)/2
                if i < d
                    % Here the inverse v-functions for the unconditional
                    % pair-copulas as part of the first tree of the C-Vine are
                    % obtained.
                    switch families(i)
                        case {3,4,5,6,12,16,17,19}
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+1));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1),rotation(i));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+1),rotation(i));
                            end
                            j = j+2;
                            
                        case {18}
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+2));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2),rotation(i));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+2),rotation(i));
                            end
                            j = j+3;
                            
                        case {0}
                            InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,0);
                            hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,0);
                            
                        otherwise
                            if rotation(i) == 0
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j));
                            else
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j),rotation(i));
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j),rotation(i));
                            end
                            j = j+1;
                            
                    end
                else
                    % Here the inverse v-functions for the second and higher
                    % trees are obtained, where for each pair-copula it is
                    % distinguished between the pair-copulas, which are partial
                    % copulas and the pair-copulas which are conditional
                    % copulas.
                    if simplified(i-d+1)
                        switch families(i)
                            case {3,4,5,6,12,16,17,19}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+1));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+1),rotation(i));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+1),rotation(i));
                                end
                                j = j+2;
                                
                            case {18}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+2));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j:j+2),rotation(i));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j:j+2),rotation(i));
                                end
                                j = j+3;
                                
                            case {0}
                                InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,0);
                                hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,0);
                                
                            otherwise
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j));
                                else
                                    InvvFunctionals{i} = @(u1,u2) VineCopulaMatlab(7,families(i),u1,u2,thetas(j),rotation(i));
                                    hFunctionals{i} = @(u1,u2) VineCopulaMatlab(4,families(i),u1,u2,thetas(j),rotation(i));
                                end
                                j = j+1;
                                
                        end
                    else
                        switch families(i)
                            case {3,4,5,6,12,16,17,19}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)]);
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)]);
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)],rotation(i));
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3)],rotation(i));
                                end
                                k = k+2;
                                
                            case {18}
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)]);
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)]);
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)],rotation(i));
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,[condparameterfunctionals{k}(u3),condparameterfunctionals{k+1}(u3),condparameterfunctionals{k+2}(u3)],rotation(i));
                                end
                                k = k+3;
                                
                            case {0}
                                InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,0);
                                hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,0);
                                
                            otherwise
                                if rotation(i) == 0
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,condparameterfunctionals{k}(u3));
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,condparameterfunctionals{k}(u3));
                                else
                                    InvvFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(7,families(i),u1,u2,condparameterfunctionals{k}(u3),rotation(i));
                                    hFunctionals{i} = @(u1,u2,u3) VineCopulaMatlab(4,families(i),u1,u2,condparameterfunctionals{k}(u3),rotation(i));
                                end
                                k = k+1;
                                
                        end
                    end
                end
            end
            
            % In the indices vector it is stored in which order the inverse
            % h-functions corresponding to the different pair-copulas have to
            % be applied to obtain a pseudo-random sample from the specified
            % C-Vine.
            indInvvfuns = zeros(1,d*(d-1)/2);
            indInvvfuns(1:3) = [1 d 2];
            indHfuns = zeros(1,d*(d-1)/2);
            indHfuns(1:3) = [1 2 d];
            
            if d > 3
                j = 4;
                for i = 3:d-1
                    indInvvfuns(j:j+i-1) = [indInvvfuns(j-i+1)+d-i+1,indInvvfuns(j-i+1:j-1)+1];
                    indHfuns(j:j+i-1) = [indHfuns(j-i+1:j-1)+1,indHfuns(j-1)+d-i+1];
                    j = j+i;
                end
            end
            
            % Here the pseudo-random sample from the specified D-Vine is
            % obtained by applying the iterative inversion algorithm for
            % generating pseudo-random variables.
            v = 1;
            h = 1;
            W = RandUniform(N,d);
            A = zeros(N,d,d);
            B = zeros(N,d,d);
            U(:,1) = W(:,1);
            A(:,1,1) = W(:,1);
            B(:,1,1) = W(:,1);
            for i = 2:d
                A(:,i,1) = W(:,i);
                for j = 2:i
                    if indInvvfuns(v)<d || simplified(indInvvfuns(v)-d+1)==1
                        A(:,i,j) = InvvFunctionals{indInvvfuns(v)}(B(:,i-1,j-1),A(:,i,j-1));
                    else
                        A(:,i,j) = InvvFunctionals{indInvvfuns(v)}(B(:,i-1,j-1),A(:,i,j-1),U(:,j:i-1));
                    end
                    v = v+1;
                end
                U(:,i) = A(:,i,i);
                B(:,i,i) = A(:,i,i);
                for j = i-1:-1:1
                    if indHfuns(h)<d || simplified(indHfuns(h)-d+1)==1
                        B(:,i,j) =  hFunctionals{indHfuns(h)}(B(:,i-1,j),A(:,i,j+1));
                    else
                        B(:,i,j) =  hFunctionals{indHfuns(h)}(B(:,i-1,j),A(:,i,j+1),U(:,j+1:i-1));
                    end
                    h = h+1;
                end
            end
    end
end

end
