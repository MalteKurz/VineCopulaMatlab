function [varargout] = VineCopulaMatlabSeedState(varargin)
%VINECOPULAMATLABSEEDSTATE Summary of this function goes here
%   Detailed explanation goes here
mlock

persistent SeedState;

if nargin==1 % The case where the state should be updated
    % check whether a valid vector is provided, i.e., a vector of 32-bit
    % unsigned integers of length 624. Normally the vector is nevertheless
    % (as usual in matlab) provided as datatype double. Meaning that Matlab
    % will check for a correct vector by transforming it into a uint32 and
    % checking whether something changed due to the transformation.
    if length(varargin{1}) ~= 624 || any(uint32(varargin{1}) ~= varargin{1})
        error('Please provide a seed state being a 32-bit unsigned integer vector of length 624. Such a seed can be obtained using the function VineCopulaMatlabSetSeed')
    end
    
    SeedState = varargin{1};
    varargout{1} = SeedState;
elseif nargin == 0 && isempty(SeedState)
    warning('VineCopulaMatlab:SetSeed','The function VineCopulaMatlabSetSeed(''shuffle'') was called to intialize a seed for the current VineCopulaMatlab session')
    SeedState = VineCopulaMatlabSetSeed('shuffle');
    varargout{1} = SeedState;
else
    varargout{1} = SeedState;
end



end

