classdef VineCopula < matlab.mixin.CustomDisplay
    %VINECOPULA The class of vine copula models (yet C-Vines and D-Vines only)
    %
    %VINECOPULA Properties
    %       dimension                - The dimension (=d) of the vine copula.
    %       type                     - The type of the vine copula, i.e.,
    %                                  'C-Vine', 'CVine' or 0 and 'D-Vine'
    %                                  ,'DVine' or 1.
    %       simplified               - A logical or vector of logicals, which
    %                                  provides the information, whether
    %                                  the simplifying assumption is
    %                                  fulfilled or not. Either the whole
    %                                  vine is set to simplified, i.e., all
    %                                  conditional copulas are chosen to be
    %                                  unconditional bivariate copulas (pair-
    %                                  copulas) or each conditional copula
    %                                  being part of the vine copula
    %                                  model is set to 1 for a pair-copula,
    %                                  i.e., a unconditional bivariate
    %                                  copula or 0 for a conditional
    %                                  copula.
    %       structure                - A vector containing the information
    %                                  about the structure of the vine
    %                                  copula model.
    %       families                 - A 2 x (d*(d-1)/2) cell array or
    %                                  matrix, where in the first row the
    %                                  copula family for each pair-copula
    %                                  is stored and in the second row, one
    %                                  can optionally store the degree of
    %                                  rotation for the pair-copula.
    %       parameters               - A 1 x (d*(d-1)/2) cell array
    %                                  containing all the parameters of the
    %                                  unconditional pair-copulas.
    %                                  Alternatively,  the same information
    %                                  (the parameters) can also be given
    %                                  simply in vectorized format, where
    %                                  the parameters are contained in the
    %                                  vector in the same ordering as for
    %                                  example the pair-copula families.
    %       condparameterfunctionals - A 1 x ((d-1)*(d-2)/2) cell array
    %                                  containing parameter functionals for
    %                                  the conditional copulas being part
    %                                  of the vine. The parameter
    %                                  functionals can also be provided in
    %                                  vectorized form.
    %       MaxLLs                   - The values of the maximized log-
    %                                  likelihoods of the PCC evaluated at
    %                                  the sequentially and jointly
    %                                  estimated parameter vectors.
    %       SeqEstParameters         - The sequentially estimated
    %                                  parameters for the vine copulas.
    %
    %VINECOPULA Methods
    %      VineCopula                - The "function" VineCopula can be
    %                                  used to construct members of the
    %                                  VineCopula class. At least the
    %                                  properties dimension, type and
    %                                  simplified have to be specified.
    %      Fit                       - The method Fit can be used to
    %                                  estimate a vine copula model,
    %                                  specified as a object of the
    %                                  VineCopula class, using a dataset by
    %                                  maximizing the log-likelihood of the
    %                                  vine copula.
    %      Sim                       - The method Sim can be used to
    %                                  simulate from a vine copula,
    %                                  specified as a object of the
    %                                  VineCopula class.
    %      StructureSelect           - The method StructureSelect can be
    %                                  used to find a "adequat" structure
    %                                  and pair-copula families (and
    %                                  parameters) for a given data set.
    %      GetPseudoObsFromVine      - The method GetPseudoObsFromVine can
    %                                  be used to obtain pseudo-
    %                                  observations from the conditional
    %                                  copulas being part of the specified
    %                                  PCC.
    %      SeqTestOnSimplified       - The method SeqTestOnSimplified can
    %                                  be used to sequentially test on the
    %                                  simplifying assumption by applying
    %                                  a vectorial independence test.
    %
    %References
    % [1]  Aas, K., C. Czado, A. Frigessi and H. Bakken (2009), "Pair-
    %      copula constructions of multiple dependence", Insurance:
    %      Mathematics and Economics 44(2), pp. 182-198.
    % [2]  Acar, E. F., C. Genest and J. Neslehová (2012), "Beyond
    %      simplified pair-copula constructions", Journal of Multivariate
    %      Analysis 110, pp. 74-90.
    % [3]  Bedford, T. and R. M. Cooke (2001), "Probability density
    %      decomposition for conditionally dependent random variables
    %      modeled by vines", Annals of Mathematics and Artificial
    %      Intelligence 32 (1), pp. 245-268.
    % [4]  Bedford, T. and R. M. Cooke (2002), "Vines -- A new graphical
    %      model for dependent random variables", The Annals of Statistics
    %      30(4), pp. 1031-1068.
    % [5]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling Dependence
    %      with C- and D-Vine Copulas: The R-Package CDVine", Journal of
    %      Statistical Software 52(3), R package version 1.1-13, pp. 1-27,
    %      url: http://CRAN.R-project.org/package=CDVine.
    % [6]  Hobaek-Haff, I., K. Aas and A. Frigessi (2010), "On the
    %      simplified pair-copula construction -- Simply useful or too
    %      simplistic?", Journal of Multivariate Analysis 101(5), pp. 1296-
    %      1310.
    % [7]  Joe, H. (1996), "Families of m-Variate Distributions With Given
    %      Margins and m(m-1)/2 Bivariate Dependence Parameters",
    %      Distributions with Fixed Marginals and Related Topics, ed. by L.
    %      Rüschendorf, B. Schweizer, and M. D. Taylor, Hayward, CA:
    %      Institute of Mathematical Statistics.
    % [8]  Joe, H. (1997), Multivariate models and dependence concepts,
    %      1. ed., reprint., Monographs on statistics and applied
    %      probability; 73, Boca Raton, Fla. [u.a.]: Chapman & Hall/CRC.
    % [9]  Kurowicka, D. and H. Joe (Eds.) (2011), "Dependece Modeling --
    %      Vine Copula Handbook", Singapore: World Scientific Publishing Co.
    %      Pte. Ltd.
    % [10] Kurz, M. (2012), "On the simplified pair-copula construction --
    %      Graphical and computational illustrations", Unpublished Term
    %      Paper, Master Seminar on Financial Econometrics, Department of
    %      Statistics, Ludwig-Maximilians-University Munich.
    % [11] Kurz, M. (2013), "Tests on the partial copula", Unpublished
    %      Master's Thesis, Department of Statistics, Ludwig-Maximilians-
    %      University Munich.
    % [12] Nikoloulopoulos, A. K., H. Joe, H. Li (2012), "Vine
    %      copulas with asymmetric tail dependence and applications
    %      to financial return data", Computational Statistics &
    %      Data Analysis 56(11), pp. 3659-3673.
    % [13] Schepsmeier, U., J. Stöber and E. C. Brechmann (2013),
    %      VineCopula: Statistical inference of vine copulas, R package
    %      version 1.2, url: http://CRAN.R-project.org/package=VineCopula.
    % [14] Spanhel, F. and M. Kurz (2014), "Simplified vine copula
    %      approximations -- Properties and consequences", submitted for
    %      publication. 
    % [15] Stöber, J., H. Joe and C. Czado (2013), "Simplified pair copula
    %      constructions -- Limitations and extensions", Journal of
    %      Multivariate Analysis 119, pp. 101-118.
    %
    %
    %
    % Author: Malte Kurz
    
    properties (SetAccess = protected)
        % dimension                  - Possible values for the dimension are
        %                              all positive integers larger than 2.
        dimension
        % type                       - Possible values for type are 'CVine',
        %                              'C-Vine', 0, 'DVine', 'D-Vine' or 1.
        type
    end
    properties
        % simplified                 - Possible values for simplified are 1
        %                              (for unconditional bivariate
        %                              copulas; pair-copulas) and 0 for
        %                              conditional copulas. If simplified
        %                              is set to 1 then every pair-copula
        %                              being part of the PCC is set to be a
        %                              partial copula. In contrast, if the
        %                              simplified property is chosen to be
        %                              0 then every conditional copula is
        %                              set to be a conditional copula with
        %                              functional parameter (note that this
        %                              can also be a unconditional copulas
        %                              whenever the functional parameter is
        %                              a constant). Alternatively, the
        %                              simplified property can also be
        %                              given as a 1 x ((d-1)*(d-2)/2)
        %                              matrix, where each entry is either
        %                              0 or 1 and therefore every copula is
        %                              set to be a pair-copula or a
        %                              conditional copula with functional
        %                              parameter, separately.
        simplified = 1
        % structure                  - Possible values for the structure in
        %                              the case of a C-Vine copula are all
        %                              permutations of the (1,...,d)
        %                              sequence. The first entry of the
        %                              structure is the root node of tree 1
        %                              (i.e., the unique node of degree d-1),
        %                              the second entry is the root node of
        %                              degree d-2 in tree 2, etc.
        %                              For a D-Vine possible values for the
        %                              structure vector are also all
        %                              permutations of the (1,...,d)
        %                              sequence. The choice [1 2 3 4 5]
        %                              corresponds to the D-Vine consisting
        %                              of the pair-copulas C_12, C_23,
        %                              C_34, C_45, C_13|2, C_24|3, C_35|4,
        %                              C_14|23, C_25|34 and C_15|234.
        structure
        % families                   - Possible values for families is either
        %                              a 1 x (d*(d-1)/2) cell array (if there
        %                              are no rotated pair-copulas within the
        %                              PCC) or a 2 x (d*(d-1)/2) cell array,
        %                              where the first row consists of the
        %                              pair-copula families and the second
        %                              one of the degrees of rotation of the
        %                              different pair-copulas. If rotation is
        %                              empty for a pair-copula, then it is
        %                              the standard (unrotated) copula of
        %                              the specified family.
        %                              The copula families and rotation
        %                              degrees can also be given in the
        %                              numerical coding (e.g. 7 for the
        %                              Clayton copula and 90 for 'r90').
        %                              The order of the entries into the
        %                              columns is tree by tree and within
        %                              the tree it is ordered according to
        %                              the structure. For example if the
        %                              structure is [2,1,3,4] then the
        %                              entries in the family cell array
        %                              correspond to the pair-copulas in the
        %                              following ordering (for the C-Vine):
        %                              C_21, C_23, C_24, C_13|2, C_14|2 and
        %                              C_34_21. Here it is important to note
        %                              that the root in each tree is always
        %                              the first variable of each pair-
        %                              copula, which is important for the
        %                              specification of pair-copulas with
        %                              asymetric dependence. For a
        %                              five-dimensional D-Vine copula with
        %                              structure [1 2 3 4 5], the ordering
        %                              of the families has to correspond to
        %                              the pair-copulas C_12, C_23, C_34,
        %                              C_45, C_13|2, C_24|3, C_35|4,
        %                              C_14|23, C_25|34 and C_15|234.
        families
        % parameters                 - Possible values for parameters are
        %                              either a 1 x (d*(d-1)/2) cell array,
        %                              where each entry is a vector of the
        %                              corresponding pair-copula
        %                              parameter(s). Alternatively, one can
        %                              also specify the parameters by
        %                              vectorizing the information contained
        %                              in the cell array of parameters.
        parameters
        % condparameterfunctionals   - Possible values for the
        %                              condparameterfunctionals are either a
        %                              1 x ((d-1)*(d-2)/2) cell array, where
        %                              each entry is a vector of the
        %                              corresponding pair-copula parameter
        %                              functional(s). Alternatively one can
        %                              also specify the parameter
        %                              functional(s) by vectorizing the
        %                              information contained in the cell
        %                              array of parameter functional(s).
        condparameterfunctionals
    end
    properties (SetAccess = protected)
        % MaxLLs                     - The first entry always corresponds to
        %                              the maximized value of the vine copula
        %                              log-likelihood evaluated at the
        %                              sequentially estimated parameters,
        %                              while the second entry corresponds to
        %                              the jointly estimated parameters.
        MaxLLs
        % SeqEstParameters           - The estimates from estimating the
        %                              vine copula by applying the
        %                              sequential estimation approach.
        SeqEstParameters
    end
    
    
    methods
        function obj = VineCopula(dimension,type,simplified,varargin)
            %VINECOPULA Constructor for objects of the VineCopula class
            % Purpose
            %        The "function" VineCopula can be used to specify
            %        members of the VineCopula class. At least the
            %        properties dimension, type and simplified have to be
            %        specified.
            %
            %
            % Usage
            %        Three entries
            %               VineCopulaObject = VineCopula(dimension,type,simplified)
            %        Four entries
            %               VineCopulaObject = VineCopula(dimension,type,simplified,structure)
            %        Five entries
            %               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families)
            %        Six entries
            %               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families,parameters)
            %        Seven entries
            %               VineCopulaObject = VineCopula(dimension,type,simplified,structure,families,parameters,condparameterfunctionals)
            %
            %        Example 1: Three-dimensional non-simplified C-Vine
            %        copula (in three different but equivalent forms)
            %               VineCopulaObject = VineCopula(3,'C-Vine',[0],[1 2 3],{'Clayton','Gumbel','Frank'},{1.2, 3,[]},{@(u) (4.*u-2).^3})
            %               VineCopulaObject = VineCopula(3,'C-Vine',[0],[1 2 3],[7, 11, 9],{1.2, 3,[]},{@(u) (4.*u-2).^3})
            %               VineCopulaObject = VineCopula(3,'C-Vine',0,[1 2 3],{'Clayton','Gumbel','Frank'},[1.2, 3],{@(u) (4.*u-2).^3})
            %        Example 2: Three-dimensional simplified C-Vine copula
            %        (in three different but equivalent forms)
            %               VineCopulaObject = VineCopula(3,'C-Vine',[1],[1 2 3],{'Clayton','Gumbel','Indep'},{1.2, 3,[]})
            %               VineCopulaObject = VineCopula(3,'C-Vine',[1],[1 2 3],[7, 11, 0],{1.2, 3,[]})
            %               VineCopulaObject = VineCopula(3,'C-Vine',1,[1 2 3],{'Clayton','Gumbel','Indep'},[1.2, 3])
            %
            % References (for the examples):
            %      Acar, E. F., C. Genest and J. Neslehová (2012), "Beyond
            %      simplified pair-copula constructions", Journal of
            %      Multivariate Analysis 110, pp. 74-90.
            %
            %
            %
            % Author: Malte Kurz
            
            narginchk(3,7)
            
            if sum(isempty(dimension)+isempty(type)+isempty(simplified)) == 0
                obj.dimension = dimension;
                obj.type = type;
                obj.simplified = simplified;
            else
                error('The properties dimension, type and simplified have to be specified.')
            end
            
            if nargin == 4
                obj.structure = varargin{1};
            elseif nargin == 5
                obj.structure = varargin{1};
                obj.families = varargin{2};
            elseif nargin == 6
                obj.structure = varargin{1};
                obj.families = varargin{2};
                obj.parameters = varargin{3};
                if sum(strcmp(obj.simplified,{'conditional'})) > 0
                    warning('VineCopula:Constructor','You have specified a non-simplified PCC, without specifying any parameter functional.')
                end
            elseif nargin == 7
                obj.structure = varargin{1};
                obj.families = varargin{2};
                obj.parameters = varargin{3};
                obj.condparameterfunctionals = varargin{4};
            elseif nargin > 7
                error('Too many input arguments.')
            end
            
        end
        
        %% Checking for the right type
        function obj = set.type(obj,type)
            if ~(strcmpi(type,'CVine') ||...
                    strcmpi(type,'C-Vine')||...
                    strcmpi(type,'DVine') ||...
                    strcmpi(type,'D-Vine') ||...
                    type == 0 ||...
                    type == 1)
                %% Not implemented yet: R-Vines
                %                      ||...
                %                     strcmpi(type,'RVine') ||...
                %                     strcmpi(type,'R-Vine'))
                error('type must be ''CVine'', ''C-Vine'', ''DVine'', ''D-Vine'', 0 (for C-Vine) or 1 (for D-Vine).')
            end
            if not(isnumeric(type))
                type = ceil(find(strcmp(type,{'CVine','C-Vine','DVine','D-Vine'}))/2)-1;
            end
            obj.type = type;
            
        end
        
        %% Checking for a allowed number of dimensions
        function obj = set.dimension(obj,dimension)
            if dimension < 3 || dimension ~= round(dimension)
                error('The dimension of the vine copula is not a natural number larger than two.')
            end
            
            obj.dimension = dimension;
            
        end
        
        
        %% Checking for the right number of indicators, which indicate,
        % whether the copulas from the second tree on are partial or
        % conditional copulas. Furthermore, checking whether all indicators
        % are either partial or conditional.
        function obj = set.simplified(obj,simplified)
            d = obj.dimension;
            
            if length(simplified)==1
                if simplified == 1
                    simplified = ones(1,(d-1)*(d-2)/2);
                elseif simplified == 0
                    simplified = zeros(1,(d-1)*(d-2)/2);
                end
            end
            if length(simplified) ~= (d-1)*(d-2)/2
                error('Wrong number of copulas are specified as unconditional or conditional.')
            elseif sum(simplified==0)+sum(simplified==1) ~= (d-1)*(d-2)/2
                error('The only possible choices for the the simplified vector are 0 (unconditional bivariate copula; pair-copula) and 1 (conditional copula).')
            end
            
            obj.simplified = simplified;
            Check(obj);
            
        end
        
        %% Check the structure of the vine
        function obj = set.structure(obj,structure)
            
            if length(structure) ~= obj.dimension
                error('Wrong number of elements in the structure vector.')
            else
                if sort(structure) ~= 1:obj.dimension
                    error('The chosen structure of the vine is not possible. Each index has to appear exactly once.')
                end
            end
            
            obj.structure = structure;
            Check(obj);
        end
        
        %% Checking for the right number of copula families
        function obj = set.families(obj,families)
            d=obj.dimension;
            
            if length(families(1,:)) ~= (d-1)*d/2
                error('Wrong number of pair-copula families.')
            end
            
            % Transfering the copula families and rotation degrees into
            % matrices.
            NumbPCs = size(families,2);
            if not(isnumeric(families))
                Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                Rotations = {'r90','r180','r270'};
                for i=1:NumbPCs
                    if sum(strcmpi(families(1,i),Families))
                        families{1,i} = find(strcmp(families(1,i),Families))-1;
                        if size(families,1) == 1 || isempty(families{2,i})
                            families{2,i} = 0;
                        else
                            families{2,i} = find(strcmp(families(2,i),Rotations)).*90;
                        end
                    else
                        error(['The copula family ' families{1,i} ' is not implemented'])
                    end
                end
                families = cell2mat(families);
            end
            
            % Checking, whether all independence copulas are specified as
            % unconditional copulas.
            if sum(min(families(1,d:end)==0,obj.simplified==0)) > 0
                error('The independence copula is always a unconditional copula and cannot be specified as a conditional copula.')
            end
            if size(families,1) > 2
                error('The families matrix has too many rows.')
            elseif size(families,1) == 1
                families = [families;zeros(size(families))];
            end
            
            % No check for combinability of family and rotation yet.
            
            obj.families = families;
            Check(obj);
        end
        
        %% Checking for the right number of parameters for the unconditional copulas
        function obj = set.parameters(obj,parameters)
            d = obj.dimension;
            Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
            
            if iscell(parameters)
                parametersVec = cell2mat(parameters);
            else
                if size(parameters,1) == 1
                    parametersVec = parameters;
                elseif size(parameters,2) == 1
                    parametersVec = parameters';
                else
                    % The case of exclusively unconditional independence
                    % copulas and in the higher order trees exclusively
                    % unconditional copulas that are either independence
                    % copulas or conditional copulas.
                    parametersVec = parameters;
                end
            end
            
            if not(isempty(parametersVec))
                
                % Checking for the right number of parameters for the
                % unconditional copulas.
                NumbParams = d*(d-1)/2 + sum(obj.families(1,:)==3) + sum(obj.families(1,:)==4) + sum(obj.families(1,:)==5) + sum(obj.families(1,:)==6) + sum(obj.families(1,:)==12) + sum(obj.families(1,:)==16) + sum(obj.families(1,:)==17) + sum(obj.families(1,:)==19)+ 2.*sum(obj.families(1,:)==18) - sum(obj.families(1,:)==0)...
                    -sum(obj.simplified==0)-sum((obj.simplified==0).*((obj.families(1,d:end)==3) + (obj.families(1,d:end)==4) + (obj.families(1,d:end)==5) + (obj.families(1,d:end)==6) + (obj.families(1,d:end)==12) + (obj.families(1,d:end)==16) + (obj.families(1,d:end)==17) + (obj.families(1,d:end)==19) + 2.*(obj.families(1,d:end)==18)));
                if length(parametersVec) == 1
                    parametersVec = repmat(parametersVec,1,NumbParams);
                elseif length(parametersVec) ~= NumbParams
                    error('Wrong number of copula parameters.')
                end
                
                families = obj.families(1,logical([ones(1,d-1) obj.simplified]));
                if size(obj.families,1) == 2
                    rotation = obj.families(2,logical([ones(1,d-1) obj.simplified]));
                end
                
                % Number of unconditional copulas
                m = length(families);
                
                [lb,ub] = PairCopulaParameterBounds(families);
                
                if sum(parametersVec < lb) || sum(parametersVec > ub)
                    k = 1;
                    for i = 1:m
                        if not(families(1,i)==0)
                            [lb,ub] = PairCopulaParameterBounds(families(1,i));
                            for j = 1:length(lb)
                                if parametersVec(k) < lb(j) || parametersVec(k) > ub(j)
                                    error(['The ' num2str(j) '. parameter of the ' Families{families(1,i)+1} ' copula has to lie between ' num2str(lb(j)) ' and ' num2str(ub(j)) '.'])
                                end
                                k = k+1;
                            end
                        end
                    end
                end
                
            else
                if not(sum(obj.simplified==0) + sum(obj.families(1,:)==0) == d*(d-1)/2)
                    error('When there is no parameter specified, then all copulas have to be independence copulas or conditional copulas with functional parameters.')
                end
                
            end
            
            if not(iscell(parameters))
                parameters = cell(1,d*(d-1)/2);
                I = 1;
                for i = 1:d*(d-1)/2
                    if i < d || obj.simplified(i-d+1)
                        switch obj.families(1,i)
                            case {3,4,5,6,12,16,17,19}
                                parameters{i} = parametersVec(I:I+1);
                                I = I+2;
                                
                            case {18}
                                parameters{i} = parametersVec(I:I+2);
                                I = I+3;
                                
                            case {0}
                                
                            otherwise
                                parameters{i} = parametersVec(I);
                                I = I+1;
                        end
                    end
                end
            end
            
            obj.parameters = parameters;
            Check(obj);
            
        end
        
        %% Checking for the right number of parameter functionals for the
        % conditional copulas.
        function obj = set.condparameterfunctionals(obj,condparameterfunctionals)
            d = obj.dimension;
            
            if not(isempty(condparameterfunctionals))
                if sum(obj.simplified==0) > 0
                    NumbParams = sum(obj.simplified==0)+sum((obj.simplified==0).*((obj.families(1,d:end)==3) + (obj.families(1,d:end)==4) + (obj.families(1,d:end)==5) + (obj.families(1,d:end)==6) + (obj.families(1,d:end)==12) + (obj.families(1,d:end)==16) + (obj.families(1,d:end)==17) + (obj.families(1,d:end)==19)+ 2.*(obj.families(1,d:end)==18)));
                    if length(condparameterfunctionals) == 1
                        condparameterfunctionalsVec = repmat(condparameterfunctionals,1,NumbParams);
                    elseif length(condparameterfunctionals)+sum(cellfun(@iscell, condparameterfunctionals)) == (d-1)*(d-2)
                        condparameterfunctionalsVec = condparameterfunctionals(not(cellfun(@isempty, condparameterfunctionals)));
                        condparameterfunctionalsVec = cat(2,condparameterfunctionalsVec{:});
                    else
                        condparameterfunctionalsVec = condparameterfunctionals;
                    end
                    if  length(condparameterfunctionalsVec) ~= NumbParams
                        error('Wrong number of parameter functionals for conditional copulas.')
                    end
                elseif sum(obj.simplified==0) == 0
                    error('Wrong number of parameter functionals for conditional copulas.')
                end
            else
                if sum(obj.simplified==0) > 0
                    error('Wrong number of parameter functionals for conditional copulas.')
                end
            end
            
            % Get the series of triangular numbers
            tri = zeros(1,d-2);
            for i = 1:(d-2)
                tri(i) = i*(i+1)/2;
            end
            
            k = 1;
            for i = 1:(d-1)*(d-2)/2
                if obj.simplified(i)==0
                    [lb,ub] = PairCopulaParameterBounds(obj.families(1,i+d-1));
                    for j = 1:length(lb)
                        condparameters = condparameterfunctionalsVec{k}(rand(500,sum(((d-1)*(d-2)/2-j)<tri)));
                        if sum(condparameters < lb(j)) > 0 || sum(condparameters > ub(j)) > 0
                            error(['The domain of the parameter functional for the ' num2str(j) '. parameter of the conditional ' Families{obj.families(1,i+d-1)} ' copula has to be a subset of the intervall (' num2str(lb(j)) ' , ' num2str(ub(j)) ').'])
                        end
                        k = k+1;
                    end
                end
            end
            
            if not(length(condparameterfunctionals)+sum(cellfun(@iscell, condparameterfunctionals)) == (d-1)*(d-2))
                condparameterfunctionals = cell(1,(d-1)*(d-2)/2);
                I = 1;
                for i = 1:(d-1)*(d-2)/2
                    if obj.simplified(i)==0
                        switch obj.families(1,i+d-1)
                            case {3,4,5,6,12,16,17,19}
                                condparameterfunctionals{i} = condparameterfunctionalsVec(I:I+1);
                                I = I+2;
                                
                            case {18}
                                condparameterfunctionals{i} = condparameterfunctionalsVec(I:I+2);
                                I = I+3;
                                
                            case {0}
                                
                            otherwise
                                condparameterfunctionals{i} = condparameterfunctionalsVec(I);
                                I = I+1;
                                
                        end
                    end
                end
            end
            
            obj.condparameterfunctionals = condparameterfunctionals;
            Check(obj);
            
        end
        
    end
    
    methods (Access = protected)
        function displayScalarObject(obj)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(obj);
            scalarHeader = className;
            header = sprintf('%s\n',scalarHeader);
            disp(header)
            propgroup = getPropertyGroups(obj);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgroup)
        end
        
        function propgrp = getPropertyGroups(obj)
            if ~isscalar(obj)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
            else
                
                Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                
                % property groups for scalars
                gTitle1 = 'Vine copula properties:';
                
                if obj.type == 0
                    Type = 'C-Vine';
                elseif obj.type == 1
                    Type = 'D-Vine';
                end
                
                propList1 = struct('Type',Type,...
                    'Dimension',obj.dimension,...
                    'Structure',obj.structure,...
                    'simplified',obj.simplified,...
                    'MaxLLs',obj.MaxLLs);
                propgrp(1) = matlab.mixin.util.PropertyGroup(propList1,gTitle1);
                
                gTitle2 = 'Building blocks tree by tree:';
                propgrp(2) = matlab.mixin.util.PropertyGroup('',gTitle2);
                
                d = obj.dimension;
                structure = obj.structure;
                
                Copula = cell(1,(d-1)*d/2);
                
                for i=1:d-1
                    Copula{i} = 'Pair_Copula';
                end
                
                for i=2:d-1
                    K= (i-1)*d-(i-1)*i/2;
                    for j=1:d-i
                        if obj.simplified(K-d+1+j) == 1
                            Copula{K+j} = 'Pair_Copula';
                        else
                            Copula{K+j} = 'Cond_Copula';
                        end
                    end
                end
                
                if obj.type == 0
                    Copula_Info = cell(1,(d-1)*d/2);
                    for i=1:d-1
                        if obj.families(2,i) == 0
                            switch obj.families(1,i)
                                case {3,4,5,6,12,16,17,19}
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) '))'];
                                    
                                case {18}
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) ','  num2str(obj.parameters{i}(3),4) '))'];
                                    
                                case {0}
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula'];
                                    
                                otherwise
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = ' num2str(obj.parameters{i}(1),4) ')'];
                            end
                        else
                            switch obj.families(1,i)
                                case {3,4,5,6,12,16,17,19}
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) '))'];
                                    
                                case {18}
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) ','  num2str(obj.parameters{i}(3),4) '))'];
                                    
                                otherwise
                                    Copula_Info{i} = ['C_' num2str(structure(1)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = ' num2str(obj.parameters{i}(1),4) ')'];
                            end
                        end
                    end
                    
                    Cond = num2str(structure(1));
                    
                    for i=2:d-1
                        K= (i-1)*d-(i-1)*i/2;
                        for j=1:d-i
                            if obj.simplified(K-d+1+j) == 1
                                if obj.families(2,K+j) == 0
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) '))'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) ','  num2str(obj.parameters{K+j}(3),4) '))'];
                                            
                                        case {0}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = ' num2str(obj.parameters{K+j}(1),4) ')'];
                                    end
                                else
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) '))'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) ','  num2str(obj.parameters{K+j}(3),4) '))'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = ' num2str(obj.parameters{K+j}(1),4) ')'];
                                    end
                                end
                            else
                                if obj.families(2,K+j) == 0
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ')'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ', ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{3}) ')'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ')'];
                                            
                                    end
                                else
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ')'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ', ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{3}) ')'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(i)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ')'];
                                            
                                    end
                                end
                            end
                        end
                        Cond = [Cond ',' num2str(structure(i))];
                    end
                    
                else
                    Copula_Info = cell(1,(d-1)*d/2);
                    for i=1:d-1
                        if obj.families(2,i) == 0
                            switch obj.families(1,i)
                                case {3,4,5,6,12,16,17,19}
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) '))'];
                                    
                                case {18}
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) ','  num2str(obj.parameters{i}(3),4) '))'];
                                    
                                case {0}
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula'];
                                    
                                otherwise
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (theta = ' num2str(obj.parameters{i}(1),4) ')'];
                            end
                        else
                            switch obj.families(1,i)
                                case {3,4,5,6,12,16,17,19}
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) '))'];
                                    
                                case {18}
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = (' num2str(obj.parameters{i}(1),4) ',' num2str(obj.parameters{i}(2),4) ','  num2str(obj.parameters{i}(3),4) '))'];
                                    
                                otherwise
                                    Copula_Info{i} = ['C_' num2str(structure(i)) ',' num2str(structure(i+1)) ': ' Families{obj.families(1,i)+1} ' copula (' num2str(obj.families(2,i)) '°, theta = ' num2str(obj.parameters{i}(1),4) ')'];
                            end
                        end
                    end
                    
                    for i=2:d-1
                        K= (i-1)*d-(i-1)*i/2;
                        for j=1:d-i
                            Cond = num2str(structure(j+1));
                            for k=j+2:i+j-1
                                Cond = [Cond ',' num2str(structure(k))];
                            end
                            if obj.simplified(K-d+1+j) == 1
                                if obj.families(2,K+j) == 0
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) '))'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) ','  num2str(obj.parameters{K+j}(3),4) '))'];
                                            
                                        case {0}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (theta = ' num2str(obj.parameters{K+j}(1),4) ')'];
                                    end
                                else
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) '))'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = (' num2str(obj.parameters{K+j}(1),4) ',' num2str(obj.parameters{K+j}(2),4) ','  num2str(obj.parameters{K+j}(3),4) '))'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) ';' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, theta = ' num2str(obj.parameters{K+j}(1),4) ')'];
                                    end
                                end
                            else
                                if obj.families(2,K+j) == 0
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ')'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ', ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{3}) ')'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ')'];
                                            
                                    end
                                else
                                    switch obj.families(1,K+j)
                                        case {3,4,5,6,12,16,17,19}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ')'];
                                            
                                        case {18}
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ', ' func2str(obj.condparameterfunctionals{K-d+1+j}{2}) ' and ' func2str(obj.condparameterfunctionals{K-d+1+j}{3}) ')'];
                                            
                                        otherwise
                                            Copula_Info{K+j} = ['C_' num2str(structure(j)) ',' num2str(structure(i+j)) '|' Cond  ': ' Families{obj.families(1,K+j)+1} ' copula (' num2str(obj.families(2,K+j)) '°, ' func2str(obj.condparameterfunctionals{K-d+1+j}{1}) ')'];
                                            
                                    end
                                end
                            end
                        end
                    end
                    
                end
                
                for i=1:d-1
                    K= (i-1)*d-(i-1)*i/2;
                    propgrp(K+3) = matlab.mixin.util.PropertyGroup(struct(Copula{K+1},Copula_Info{K+1}),[num2str(i) '. Tree']);
                    for j=2:d-i
                        propgrp(K+j+2) = matlab.mixin.util.PropertyGroup(struct(Copula{K+j},Copula_Info{K+j}));
                    end
                end
                
            end
        end
    end
    
    methods (Access = private)
        function check = Check(obj)
            %CHECK Check objects of the VineCopula class for correctly specified properties
            % Purpose
            %        The private function is designed for cross-checking
            %        the specified vine copula properties. Therefore, it is
            %        called whenever any property is set or changed. There
            %        is no need to call this function being a user of the
            %        toolbox, as it is always automatically called whenever
            %        a object of the class VineCopula is generated or
            %        properties are changed.
            %
            %
            % Usage
            %        Checking the properties of a VineCopula object
            %               check = Check(VineCopula)
            %
            %
            % Inputs
            %        VineCopula   = An object from the class VineCopula.
            %
            %
            % Outputs
            %       check         = A logical, which is 1/true whenver the
            %                       properties are set correctly.
            %
            %
            %
            % Author: Malte Kurz
            
            d = obj.dimension;
            Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
            
            % Checking for the right number of indicators, which indicate,
            % whether the copulas from the second tree on are partial or
            % conditional copulas. Furthermore, checking whether all indicators
            % are either partial or conditional.
            if length(obj.simplified) ~= (d-1)*(d-2)/2
                warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of copulas are specified as partial or conditional.')
            elseif sum(obj.simplified==0)+sum(obj.simplified==1) ~= (d-1)*(d-2)/2
                warning('VineCopula:ChangeOfOtherPropertyNeeded','The only possible choices for the the simplified vector are 0 (unconditional bivariate copula; pair-copula) and 1 (conditional copula).')
            end
            
            % Check the structure of the vine
            if not(isempty(obj.structure))
                if length(obj.structure) ~= obj.dimension
                    warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of elements in the structure vector.')
                else
                    if sort(obj.structure) ~= 1:obj.dimension
                        warning('VineCopula:ChangeOfOtherPropertyNeeded','The chosen structure of the vine is not possible.')
                    end
                end
            end
            
            % Checking for the right number of copula families
            if not(isempty(obj.families))
                if length(obj.families(1,:)) ~= (d-1)*d/2
                    warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of pair-copula families.')
                end
                % Checking, whether all independence copulas are specified as
                % partial copulas.
                if sum(min(obj.families(1,d:end)==0,obj.simplified==0)) > 1
                    warning('VineCopula:ChangeOfOtherPropertyNeeded','The independence copula is always a unconditional copula and cannot be specified as a conditional copula.')
                end
                if size(obj.families,1) > 2
                    error('The families cell array has too many rows.')
                end
                
            end
            
            % Checking for the right number of parameters for the unconditional copulas
            if not(isempty(obj.parameters))
                
                parameters = cell2mat(obj.parameters);
                
                % Checking for the right number of parameters for the partial and
                % unconditional copulas.
                if length(parameters) ~= d*(d-1)/2 + sum(obj.families(1,:)==3) + sum(obj.families(1,:)==4) + sum(obj.families(1,:)==5) + sum(obj.families(1,:)==6) + sum(obj.families(1,:)==12) + sum(obj.families(1,:)==16) + sum(obj.families(1,:)==17) + sum(obj.families(1,:)==19)+ 2.*sum(obj.families(1,:)==18) - sum(obj.families(1,:)==0)...
                        -sum(obj.simplified==0)-sum((obj.simplified==0).*((obj.families(1,d:end)==3) + (obj.families(1,d:end)==4) + (obj.families(1,d:end)==5) + (obj.families(1,d:end)==6) + (obj.families(1,d:end)==12) + (obj.families(1,d:end)==16) + (obj.families(1,d:end)==17) + (obj.families(1,d:end)==19)+ 2.*(obj.families(1,d:end)==18)))
                    warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of copula parameters.')
                end
                
                d = obj.dimension;
                
                families = obj.families(1,logical([ones(1,d-1) obj.simplified]));
                if size(obj.families,1) == 2
                    rotation = obj.families(2,logical([ones(1,d-1) obj.simplified]));
                end
                
                % Number of unconditional and partial copulas
                m = length(families);
                
                [lb,ub] = PairCopulaParameterBounds(families);
                
                if isempty(parameters)
                    if not(length(families(1,:)) == sum(families(1,:)==0))
                        warning('VineCopula:ChangeOfOtherPropertyNeeded','When there is no parameter specified, then all copulas have to be independence copulas or conditional copulas with functional parameter.')
                    end
                else
                    if sum(parameters < lb) || sum(parameters > ub)
                        k = 1;
                        for i = 1:m
                            if not(families(1,i)==0)
                                [lb,ub] = PairCopulaParameterBounds(families(1,i));
                                for j = 1:length(lb)
                                    if parameters(k) < lb(j) || parameters(k) > ub(j)
                                        warning('VineCopula:ChangeOfOtherPropertyNeeded',['The ' num2str(j) '. parameter of the ' Families{obj.families(1,i)} ' copula has to lie between ' num2str(lb(j)) ' and ' num2str(ub(j)) '.'])
                                    end
                                    k = k+1;
                                end
                            end
                        end
                    end
                end
                
            end
            
            % Checking for the right number of parameter functionals for the
            % conditional copulas.
            if not(isempty(obj.condparameterfunctionals))
                condparameterfunctionals = obj.condparameterfunctionals(not(cellfun(@isempty, obj.condparameterfunctionals)));
                condparameterfunctionals = cat(2,condparameterfunctionals{:});
                if sum(obj.simplified==0) > 0
                    if length(condparameterfunctionals) ~= sum(obj.simplified==0)+sum((obj.simplified==0).*((obj.families(1,d:end)==3) + (obj.families(1,d:end)==4) + (obj.families(1,d:end)==5) + (obj.families(1,d:end)==6) + (obj.families(1,d:end)==12) + (obj.families(1,d:end)==16) + (obj.families(1,d:end)==17) + (obj.families(1,d:end)==19)+ 2.*(obj.families(1,d:end)==18)))
                        warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of parameter functionals for conditional copulas.')
                    end
                elseif sum(strcmp(obj.simplified,{'conditional'})) == 0
                    warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of parameter functionals for conditional copulas.')
                end
                % Hided, because warnings are coming if one specifies a
                % vine with a conditional parameter functional.
                %else
                %    if sum(obj.simplified==0)>0
                %        warning('VineCopula:ChangeOfOtherPropertyNeeded','Wrong number of parameter functionals for conditional copulas')
                %    end
                
                % Get the series of triangular numbers
                tri = zeros(1,d-2);
                for i = 1:(d-2)
                    tri(i) = i*(i+1)/2;
                end
                
                k = 1;
                for i = 1:(d-1)*(d-2)/2
                    if obj.simplified(i)==0
                        [lb,ub] = PairCopulaParameterBounds(obj.families(1,i+d-1));
                        for j = 1:length(lb)
                            condparameters = condparameterfunctionals{k}(rand(500,sum(((d-1)*(d-2)/2-j)<tri)));
                            if sum(condparameters < lb(j)) > 0 || sum(condparameters > ub(j)) > 0
                                error(['The domain of the parameter functional for the ' num2str(j) '. parameter of the conditional ' Families(obj.families(1,i+d-1)+1) ' copula has to be a subset of the intervall (' num2str(lb(j)) ' , ' num2str(ub(j)) ').'])
                            end
                            k = k+1;
                        end
                    end
                end
            end
            
            check = true;
            
        end
        
    end
    
    methods
        function VineCopulaHat = Fit(obj,u,varargin)
            %FIT Estimating objects of the VineCopula class
            % Purpose
            %        The function computes ML-estimates for the parameters
            %        of a simplified vine copula. Therefore, first starting
            %        values for the joint estimation are obtained by
            %        iteratively estimating the pair-copulas in the first
            %        trees and using those estimates to obtain the
            %        arguments for the copulas in the second tree. Then the
            %        pair-copulas in the second tree are estimated and so
            %        on. These estimated parameters from the sequential
            %        procedure are then used to obtained the ML-estimates,
            %        by minimizing the overall negative log-likelihood of
            %        the whole vine copula numerically.
            %
            %
            % Usage
            %        Estimating a simplified vine copula (joint estimation;
            %        the default method)
            %          VineCopulaHat = Fit(VineCopulaObject,u)
            %          VineCopulaHat = Fit(VineCopulaObject,u,'joint')
            %        Estimating a simplified vine copula (sequential
            %        estimation)
            %          VineCopulaHat = Fit(VineCopulaObject,u,'sequential')
            %        Estimating a simplified vine copula (with a cut off
            %        tree / truncation level)
            %          VineCopulaHat = Fit(VineCopulaObject,u,EstMethod,CutOffTree)
            %
            %
            % Inputs
            %        VineCopulaObject= An object from the class VineCopula.
            %        u               = A (n x d) dimensional vector of
            %                          values lying in [0,1] (the
            %                          observations).
            %        EstMethod       = The estimation method must be either
            %                          'joint' or 'sequential'. If it is
            %                          not explicitly given, a joint
            %                          estimation is performed (default).
            %        CutOffTree      = The CutOffTree (or also called
            %                          truncation level) can be used to set
            %                          all pair-copulas from the (CutOffTree
            %                          + 1)-th tree on to independence
            %                          copulas (i.e., ignore them in the
            %                          joint estimation). The CutOffTree
            %                          does only influence the joint
            %                          estimation.
            %
            %
            % Outputs
            %       VineCopulaHat   = An object from the class VineCopula.
            %                         The sequential estimates are stored
            %                         in VineCopulaHat.SeqEstParameters,
            %                         the estimated parameters from the
            %                         joint estimation are stored in
            %                         VineCopulaHat.parameters and the two
            %                         maximized values of the vine copula
            %                         log-likelihood are stored in
            %                         VineCopulaHat.MaxLLs.
            %
            %
            %
            % Author: Malte Kurz
            
            % Check the (Copula-)data input.
            CheckData(u)
            
            narginchk(2,4)
            
            if nargin == 3
                if isempty(varargin{1})
                    EstMethod = 'joint';
                else
                    EstMethod = varargin{1};
                end
            elseif nargin == 4
                if isempty(varargin{1})
                    EstMethod = 'joint';
                else
                    EstMethod = varargin{1};
                end
                if isempty(varargin{2})
                    % Set CutOffTree to dimension-1 to obtain the same
                    % results like without specifying any CutOffTree
                    CutOffTree = obj.dimension -1;
                else
                    CutOffTree = varargin{2};
                end
            else
                EstMethod = 'joint';
            end
            
            u = u(:,obj.structure);
            if size(obj.families,1) == 1
                fams = obj.families;
                if nargin == 4
                    [ParamHat, MaxLogLikes, theta0, newfams] = VineCopulaFit(obj.type,fams,obj.dimension,u,[],EstMethod,CutOffTree);
                else
                    [ParamHat, MaxLogLikes, theta0] = VineCopulaFit(obj.type,fams,obj.dimension,u,[],EstMethod);
                end
            elseif size(obj.families,1) == 2
                fams = obj.families(1,1:end);
                rotation = obj.families(2,:);
                if nargin == 4
                    [ParamHat, MaxLogLikes, theta0, newfams] = VineCopulaFit(obj.type,fams,obj.dimension,u,rotation,EstMethod,CutOffTree);
                else
                    [ParamHat, MaxLogLikes, theta0] = VineCopulaFit(obj.type,fams,obj.dimension,u,rotation,EstMethod);
                end
            end
            
            if sum(obj.simplified==0) > 0
                warning('VineCopula:Fit','The vine copula has been estimated as a simplified vine copula.')
            end
            
            if nargin == 4
                VineCopulaHat = VineCopula(obj.dimension,obj.type,1,obj.structure,newfams,ParamHat);
            else
                VineCopulaHat = VineCopula(obj.dimension,obj.type,1,obj.structure,obj.families,ParamHat);
                newfams = obj.families;
            end
            
            VineCopulaHat.MaxLLs = MaxLogLikes;
            
            d = obj.dimension;
            SeqParameters = cell(1,d*(d-1)/2);
            I = 1;
            for i = 1:d*(d-1)/2
                switch newfams(1,i)
                    case {3,4,5,6,12,16,17,19}
                        SeqParameters{i} = theta0(I:I+1);
                        I = I+2;
                        
                    case {18}
                        SeqParameters{i} = theta0(I:I+2);
                        I = I+3;
                        
                    case {0}
                        
                    otherwise
                        SeqParameters{i} = theta0(I);
                        I = I+1;
                        
                end
            end
            
            VineCopulaHat.SeqEstParameters = SeqParameters;
            
        end
        
        function U = Sim(obj,N)
            %SIM Generating (pseudo-)random observations from a vine copula
            % Purpose
            %        The function draws N pseudo-random d-dimensional
            %        tuples from a vine copula. The function can be used to
            %        sample from arbitrarily dimensional C- and D-Vines of
            %        the simplified and non-simplified form. In the case of
            %        a non-simplified C- or D-Vine, the parameters of the
            %        conditional copulas have to be given as functionals of
            %        the variables on which the conditioning is done.
            %
            %
            % Usage
            %       Simulate N observations from a vine copula
            %                      U = Sim(VineCopulaObject,N)
            %
            %
            % Inputs
            %        VineCopulaObject= An object from the class VineCopula.
            %        N               = The number of observations that
            %                          should be simulated from the vine
            %                          copula.
            %
            %
            % Outputs
            %        U               = A (N x d) Matrix of simulated tuples
            %                          from the specified vine copula, where
            %                          every row is a simulated d-dimensional
            %                          tuple from the d-dimensional vine
            %                          copula.
            %
            %
            %
            % Author: Malte Kurz
            
            parameters = cell2mat(obj.parameters);
            
            if size(obj.families,1) == 1
                if sum(obj.simplified==0) == 0
                    U = VineCopulaRand(obj.type,N,obj.dimension,obj.families,parameters,obj.simplified);
                else
                    condparameterfunctionals = obj.condparameterfunctionals(not(cellfun(@isempty, obj.condparameterfunctionals)));
                    condparameterfunctionals = cat(2,condparameterfunctionals{:});
                    U = VineCopulaRand(obj.type,N,obj.dimension,obj.families,parameters,obj.simplified,condparameterfunctionals);
                end
            elseif size(obj.families,1) == 2
                if sum(obj.simplified==0) == 0
                    U = VineCopulaRand(obj.type,N,obj.dimension,obj.families(1,:),parameters,obj.simplified,[],obj.families(2,:));
                else
                    condparameterfunctionals = obj.condparameterfunctionals(not(cellfun(@isempty, obj.condparameterfunctionals)));
                    condparameterfunctionals = cat(2,condparameterfunctionals{:});
                    U = VineCopulaRand(obj.type,N,obj.dimension,obj.families(1,:),parameters,obj.simplified,condparameterfunctionals,obj.families(2,:));
                end
            end
            
        end
        
        function NewObj = StructureSelect(obj,data,varargin)
            %STRUCTURESELECT Selecting the structure and pair-copula families for a vine copula
            % Purpose
            %        The function can be used to find an adequate structure
            %        and to select pair-copula families for a given data
            %        set. By default (StructuringRule = 0) the nodes of
            %        the C-vine are chosen in a way that in each tree,
            %        the root (i.e. the node, which is connected by a
            %        copula to all other nodes) is the variable, which has
            %        maximal dependence with all other variables. The
            %        maximal dependence is found by choosing the variable,
            %        which has the maximal column sum in the matrix of
            %        absolute empirical Kendall’s τ (cf. Schepsmeier,
            %        Stöber, and Brechmann (2013) for an R-function
            %        (RVineStructureSelect) of exactly the same procedure
            %        and Czado, Schepsmeier, and Min (2012, p. 240) for the
            %        theoretical background of the approach).
            %        Alternative methods for choosing a structure for
            %        C-Vine copulas (cf. Nikoloulopoulos, Joe and Li
            %        (2012)):
            %        The root node of the first tree is chosen in same way
            %        as in the default method, i.e., the variable with the
            %        strongest dependence with all other variables. Then
            %        one can choose between three rules suggested in
            %        Nikoloulopoulos, Joe and Li (2012, p.3665):
            %        StructuringRule = 1: List the other variables by their
            %                             dependence to the root node in
            %                             decreasing order.
            %        StructuringRule = 2: List the other variables by their
            %                             dependence to the root node in
            %                             increasing order.
            %        StructuringRule = 3: List the other variables
            %                             sequentially by choosing the
            %                             variable which is least dependent
            %                             with the previously selected one.
            %        Furthermore, the copula families are chosen
            %        according to the AIC criterion and for each pair-
            %        copula an independence test is performed (cf.
            %        Schepsmeier, Stöber, and Brechmann (2013) and
            %        Brechmann and Schepsmeier (2013) for R-functions
            %        (RVineStructureSelect / RVineCopSelect /
            %        CDVineCopSelect) and Genest and Favre (2007, p. 351)
            %        for the independence test).
            %        For D-Vine copulas there is no structuring rule
            %        implemented yet. Therefore, the function uses the
            %        specified structure and only selects the pair-copula
            %        families using the AIC criterion.
            %
            %
            % Usage
            %        Select a simplified vine copula model
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u)
            %        Select a simplified vine copula model, where the
            %        pair-copulas are from a specified set of families
            %        (possible choices 'all' (default), 'R', 'R-package',
            %        'VineCopulaPackage' (they all are equivalent and
            %        correspond to the set of pair-copulas of the R-package
            %        VineCopula (Schepsmeier, Stöber, and Brechmann (2013),
            %        Version 1.2)) or a cell-array consisting of possible
            %        pair-copula families (i.e., a user selected list of
            %        possible pair-copula families).
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,'all')
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,'R')
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,'R-package')
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,'VineCopulaPackage')
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset)
            %        Select a simplified vine copula model using an
            %        alternative structuring rule:
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,1)
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,2)
            %          VineCopulaHat = StructureSelect(VineCopulaObject,u,familyset,3)
            %
            %
            %
            % Inputs
            %        VineCopulaObject= An object from the class VineCopula.
            %        u               = A (n x d) dimensional vector of
            %                          values lying in [0,1] (the
            %                          observations).
            %        familyset       = The set of possible pair-copula
            %                          families. By setting it to the
            %                          strings 'all', 'R', 'R-package' or
            %                          'VineCopulaPackage' one can choose
            %                          one of the pre-defined sets.
            %                          Alternatively one can choose an
            %                          array containing a subset of the
            %                          possible families:
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
            % Outputs
            %        VineCopulaHat   = An object from the class VineCopula.
            %                          The select vine copula structure can
            %                          be found in VineCopulaHat.structure
            %                          and the selected pair-copulas in
            %                          VineCopulaHat.families. Furthermore,
            %                          the sequential estimates, which are
            %                          obtained during the selection
            %                          procedure of the structure and
            %                          copula families are stored in
            %                          VineCopulaHat.parameters.
            %
            %
            % References
            % [1]  Brechmann, E. C. and U. Schepsmeier (2013), "Modeling
            %      Dependence with C- and D-Vine Copulas: The R-Package
            %      CDVine", Journal of Statistical Software 52(3), R
            %      package version 1.1-13, pp. 1-27, url:
            %      http://CRAN.R-project.org/package=CDVine.
            % [2]  Czado, C., U. Schepsmeier, and A. Min (2012), "Maximum
            %      likelihood estimation of mixed C-vines with application
            %      to exchange rates", Statistical Modelling 12(3), pp.
            %      229-255.
            % [3]  Genest, C. and A. Favre (2007), "Everything You Always
            %      Wanted to Know about Copula Modeling but Were Afraid to
            %      Ask", Journal of Hydrologic Engineering 12(4), pp. 347-
            %      368.
            % [4]  Nikoloulopoulos, A. K., H. Joe, H. Li (2012), "Vine
            %      copulas with asymmetric tail dependence and applications
            %      to financial return data", Computational Statistics &
            %      Data Analysis 56(11), pp. 3659-3673.
            % [5]  Schepsmeier, U., J. Stöber, and E. C. Brechmann (2013),
            %      VineCopula: Statistical inference of vine copulas, R
            %      package version 1.2, url:
            %      http://CRAN.R-project.org/package=VineCopula.
            %
            %
            %
            % Author: Malte Kurz
            
            d = obj.dimension;
            families = zeros(2,(d-1)*d/2);
            
            % Check the (Copula-)data input.
            CheckData(data)
            
            switch obj.type
                case 0
                    narginchk(2,4)
                    
                    if nargin == 3
                        if isempty(varargin{1})
                            Familyset = -1;
                        else
                            Familyset = varargin{1};
                            if not(isnumeric(Familyset))
                                if not(iscell(Familyset))
                                    switch Familyset
                                        case {'all'}
                                            Familyset = -1;
                                            %families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                            
                                        case {'R','R-package','VineCopulaPackage'}
                                            Familyset = [3;4;5;6;7;9;10;11;13;16;17;19];
                                            %families = {'t','Gaussian','Frank','Clayton','Gumbel','Joe','BB1','BB6','BB7','BB8','Tawn1','Tawn2'};
                                            
                                    end
                                else
                                    Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                    for i=1:length(Familyset)
                                        if sum(strcmpi(Familyset(i),Families))
                                            Familyset{i} = find(strcmp(Familyset(i),Families))-1;
                                        else
                                            error(['The copula family ' Familyset{i} ' is not implemented'])
                                        end
                                    end
                                    Familyset = cell2mat(Familyset);
                                end
                            end
                        end
                        StructuringRule = 0;
                    elseif nargin == 4
                        if isempty(varargin{1})
                            Familyset = -1;
                        else
                            Familyset = varargin{1};
                            if not(isnumeric(Familyset))
                                if not(iscell(Familyset))
                                    switch Familyset
                                        case {'all'}
                                            Familyset = -1;
                                            %families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                            
                                        case {'R','R-package','VineCopulaPackage'}
                                            Familyset = [3;4;5;6;7;9;10;11;13;16;17;19];
                                            %families = {'t','Gaussian','Frank','Clayton','Gumbel','Joe','BB1','BB6','BB7','BB8','Tawn1','Tawn2'};
                                            
                                    end
                                else
                                    Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                    for i=1:length(Familyset)
                                        if sum(strcmpi(Familyset(i),Families))
                                            Familyset{i} = find(strcmp(Familyset(i),Families))-1;
                                        else
                                            error(['The copula family ' Familyset{i} ' is not implemented'])
                                        end
                                    end
                                    Familyset = cell2mat(Familyset);
                                end
                            end
                        end
                        if isempty(varargin{2})
                            StructuringRule = 0;
                        else
                            StructuringRule = varargin{2};
                            if StructuringRule ~= 0 && StructuringRule ~= 1 && StructuringRule ~= 2 && StructuringRule ~= 3
                                error('StructuringRule has to be 0, 1, 2 or 3')
                            end
                        end
                    else
                        Familyset =-1;
                        StructuringRule = 0;
                    end
                    
                    [structure,families(1,:),families(2,:),parameters]=VineStructureSelect(0,data,StructuringRule,Familyset);
                    structure = structure + 1;
                    
                    if sum(obj.simplified==0) > 0
                        warning('VineCopula:StructureSelect','The vine copula structue has been selected and estimated as a simplified vine copula.')
                    end
                    
                    NewObj = VineCopula(obj.dimension,obj.type,1,structure,families,parameters);
                    NewObj.SeqEstParameters = parameters;
                    
                case 1
                    %warning('An automatic selection of a structure is not yet implemented for D-Vines. The function uses the given structure (or by default 1:d) and only selects pair-copula families.')
                    narginchk(2,3)
                    
                    if nargin == 3
                        if isempty(varargin{1})
                            Familyset = -1;
                        else
                            Familyset = varargin{1};
                            if not(isnumeric(Familyset))
                                if not(iscell(Familyset))
                                    switch Familyset
                                        case {'all'}
                                            Familyset = -1;
                                            %families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                            
                                        case {'R','R-package','VineCopulaPackage'}
                                            Familyset = [3;4;5;6;7;9;10;11;13;16;17;19];
                                            %families = {'t','Gaussian','Frank','Clayton','Gumbel','Joe','BB1','BB6','BB7','BB8','Tawn1','Tawn2'};
                                            
                                    end
                                else
                                    Families = {'Indep','AMH','AsymFGM','BB1','BB6','BB7','BB8','Clayton','FGM','Frank','Gaussian','Gumbel','IteratedFGM','Joe','PartialFrank','Plackett','Tawn1','Tawn2','Tawn','t'};
                                    for i=1:length(Familyset)
                                        if sum(strcmpi(Familyset(i),Families))
                                            Familyset{i} = find(strcmp(Familyset(i),Families))-1;
                                        else
                                            error(['The copula family ' Familyset{i} ' is not implemented'])
                                        end
                                    end
                                    Familyset = cell2mat(Familyset);
                                end
                            end
                        end
                    else
                        Familyset =-1;
                    end
                    
                    if isempty(obj.structure)
                        obj.structure = 1:d;
                    end
                    
                    [~,families(1,:),families(2,:),parameters]=VineStructureSelect(1,data(:,obj.structure),0,Familyset);
                    
                    if sum(obj.simplified==0) > 0
                        warning('VineCopula:StructureSelect','The vine copula structue has been selected and estimated as a simplified vine copula.')
                    end
                    
                    NewObj = VineCopula(obj.dimension,obj.type,1,obj.structure,families,parameters);
                    NewObj.SeqEstParameters = parameters;
                    
            end
            
        end
        
        function [pVals,TestStats,BootTestStats] = SeqTestOnSimplified(obj,data,varargin)
            %SEQTESTONSIMPLIFIED Sequentially testing the simplified assumption for vine copulas
            % Purpose
            %        The function can be used to sequentially test the
            %        simplifying assumption for vine copulas. The procedure
            %        heavily relies on a stochastic representation of the
            %        simplifying assumption stated in Spanhel and Kurz
            %        (2014). The stochastic representation allows to test
            %        the simplifying assumption by using tests on vectorial
            %        independencies. Note, that the performed tests are 
            %        always based on the assumption that one is able to
            %        observe (pseudo-)observations from the conditional
            %        copulas that should be tested. Therefore, for
            %        interpreting the test results for a specific tree of
            %        the vine copula one needs to assume that the lower
            %        trees (including the assumption of unconditional
            %        copulas) are correctly specified.
            %
            %
            % Usage
            %        Testing the simplifying assumption sequentially
            %         [pVals,TestStats,BootTestStats] = SeqTestOnSimplified(VineCopulaObject,data,N)
            %        Testing the simplifying assumption sequentially as
            %        goodness-of-fit test (i.e., without reestimating the
            %        whole vine copula for each tree)
            %         [pVals,TestStats,BootTestStats] = SeqTestOnSimplified(VineCopulaObject,data,N,true)
            %
            %
            % Inputs
            %        VineCopulaObject      = An object from the class
            %                                VineCopula.
            %        data                  = A (n x d) dimensional vector
            %                                of values lying in [0,1] (the
            %                                observations).
            %        N                     = The number of boostrap
            %                                replications for the
            %                                multiplier bootstrap in the
            %                                vectorial independence test.
            %        GoF                   = A logical, which is by default
            %                                false. Then, before the tests
            %                                on vectorial independence are
            %                                applied in each tree, the
            %                                whole vine copula model is
            %                                reestimated up to this tree
            %                                and then the tests are
            %                                performed. Otherwise, i.e., if
            %                                it is applied as a goodness-
            %                                of-fit test, the joint
            %                                estimates given as inputs
            %                                through the VineCopula object
            %                                are used to obtain the
            %                                (pseudo-)observations in each
            %                                tree (without re-estimation).
            %
            %
            % Outputs
            %        pVal                  = A vector of p-values of the
            %                                vectorial independence tests. 
            %                                Every entry corresponds to one
            %                                specific copula being part
            %                                of the whole vine copula.
            %        TestStat              = A vector of realized values
            %                                for the test statistics.
            %                                of the whole vine copula.
            %        BootTestStats         = A matrix of realized values
            %                                for the bootstrapped test
            %                                statistics.
            %
            %
            %
            % Author: Malte Kurz
            
            % Check the (Copula-)data input.
            CheckData(data)
            
            if nargin == 4
                if not(isempty(varargin{2}))
                    GoF = varargin{2};
                else
                    GoF = false;
                end
                if not(isempty(varargin{1}))
                    N = varargin{1};
                else
                    N = 1000;
                end
            elseif nargin == 3
                GoF = false;
                if not(isempty(varargin{1}))
                    N = varargin{1};
                else
                    N = 1000;
                end
            else
                GoF = false;
                N=1000;
            end
            
            d = obj.dimension;
            
            if GoF
                % Sequentially testing the simplifying assumption as
                % goodness-of-fit test: Testing on partial copulas by using the
                % given copula parameters (i.e., in most of the cases their
                % joint estimates) to compute (pseudo-)observations.
                U = GetPseudoObsFromVine(obj,data);
                
            else
                % Sequentially testing the simplifying assumption, by
                % estimating the vine (d-2)-times in its truncated version
                % to obtain pseudo-observations. That means, that one does
                % not assume any parameteric form for the copulas that
                % should be tested against partial copulas. But, note that
                % one has to assume partial copulas for the lower trees to
                % be able to estimate the truncated vine, which is needed
                % to be able to obtain pseudo-observations from the
                % conditional copulas in the higher order trees of the
                % vine.
                
                SeqObj = obj;
                SeqObj.parameters = obj.SeqEstParameters;
                V = GetPseudoObsFromVine(SeqObj,data);
                switch obj.type
                    case 0
                        U = zeros(size(data,1),d*(d-1)/2-1);
                        U(:,1:d-1) = V(:,1:d-1);
                        F = d;
                        for j = 2:d-2
                            SeqObj = Fit(obj,data,'joint',j);
                            V = GetPseudoObsFromVine(SeqObj,data);
                            U(:,F:F+d-j-1) = V(:,F:F+d-j-1);
                            F = F+d-j;
                        end
                    case 1
                        U = zeros(size(data,1),(d-2)*(d-1)/2,2);
                        U(:,1:d-2,:) = V(:,1:d-2,:);
                        F = d;
                        for j = 2:d-2
                            SeqObj = Fit(obj,data,'joint',j);
                            V = GetPseudoObsFromVine(SeqObj,data);
                            U(:,F-j+1:F+d-2*j-1) = V(:,F-j+1:F+d-2*j-1);
                            F = F+d-j;
                        end
                end
            end
            
            data = data(:,obj.structure);
            
            pVals = ones((d-1)*(d-2)/2,1)*-1234;
            TestStats = ones((d-1)*(d-2)/2,1)*-1234;
            BootTestStats = ones(N,(d-1)*(d-2)/2);
            
            switch obj.type
                case 0
                    F = 1;
                    for j = 1:d-2
                        for i = 1:d-j-1
                            [pVals(F+i-1),TestStats(F+i-1),BootTestStats(:,F+i-1)] = VecIndepTest(U(:,[F+j-1,F+i+j-1]),data(:,1:j),N);
                        end
                        F = F+d-j-1;
                    end
                case 1
                    F = 1;
                    for j = 1:d-2
                        for i = 1:d-j-1
                            [pVals(F+i-1),TestStats(F+i-1),BootTestStats(:,F+i-1)] = VecIndepTest([U(:,F+i-1,1),U(:,F+i-1,2)],data(:,i+1:i+j),N);
                        end
                        F = F+d-j-1;
                    end
            end
            
        end
        
        function [U] = GetPseudoObsFromVine(obj,u)
            %GETPSEUDOOBSFROMVINE Computing (pseudo-)observations from vine copulas
            % Purpose
            %        The function can be used to obtain
            %        (pseudo-)observations from all pair-copulas
            %        being building blocks of a specific vine copula.
            %
            %
            % Usage
            %                      U = GetPseudoObsFromVine(VineCopulaObject,u)
            %
            %
            % Inputs
            %        VineCopulaObject= An object from the class VineCopula.
            %        u               = A (n x d) dimensional vector of
            %                          values lying in [0,1] (the
            %                          observations).
            %
            %
            % Outputs
            %        U               = A (n x (d*(d-1)/2-1)) dimensional
            %                          vector of values lying in [0,1] (the
            %                          (pseudo-)observations from the
            %                          conditional copulas being building
            %                          blocks of the vine copula). The
            %                          order of the columns is (at the
            %                          example of the 5-dimensional C-vine
            %                          with structure [1 2 3 4 5]):
            %                          2|1, 3|1, 4|1, 5|1, 3|12, 4|12,
            %                          5|12, 4|123, 5|123.
            %                          For D-Vines it is a (n x
            %                          ((d-2)*(d-1)/2) x 2) dimensional
            %                          matrix of values lying in [0,1] (the
            %                          (pseudo-)observations from the
            %                          conditional copulas being building
            %                          blocks of the vine copula). For the
            %                          5-dimensional D-Vine with structure
            %                          [1 2 3 4 5] in U(:,:,1) the order
            %                          of the columns is: 1|2, 2|3, 3|4,
            %                          1|23, 2|34, 1|234 and in U(:,:,2) it
            %                          is: 3|2, 4|3, 5|4, 4|23, 5|34, 5|234.
            %
            %
            %
            % Author: Malte Kurz
            
            % Check the (Copula-)data input.
            CheckData(u)
            
            u = u(:,obj.structure);
            thetas = cell2mat(obj.parameters);
            d = obj.dimension;
                
            if size(obj.families,1) == 1
                fams = obj.families;
                switch obj.type
                    case 0
                        U = VineGetPseudoObs(obj.type,u,fams,thetas);
                    case 1
                        U = zeros(size(u,1),(d-2)*(d-1)/2,2);
                        [U(:,:,1),U(:,:,2)] = VineGetPseudoObs(obj.type,u,fams,thetas);
                end
            elseif size(obj.families,1) == 2
                fams = obj.families(1,1:end);
                rotation = obj.families(2,:);
                switch obj.type
                    case 0
                        U = VineGetPseudoObs(obj.type,u,fams,thetas,rotation);
                    case 1
                        U = zeros(size(u,1),(d-2)*(d-1)/2,2);
                        [U(:,:,1),U(:,:,2)] = VineGetPseudoObs(obj.type,u,fams,thetas,rotation);
                end
            end
            
        end
        
    end
    
end
