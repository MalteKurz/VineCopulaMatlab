function SeedState = VineCopulaMatlabSetSeed(varargin)
%VINECOPULAMATLABSETSEED 
%   

if nargin==0 || (ischar(varargin{1}) && strcmp(varargin{1},'shuffle'))
    % Call C++ function VineCopulaMatlabSetSeed
    SeedState = VineCopulaMatlab(1005);
elseif nargin==1 && not(mod(varargin{1},1))
    SeedState = VineCopulaMatlab(1005,varargin{1});
end

VineCopulaMatlabSeedState(SeedState);


end

