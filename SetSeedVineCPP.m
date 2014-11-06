function SetSeedVineCPP
%SETSEEDVINECPP Setting the global seed for random number generation in C++
% Purpose
%        The function can be used to set the global (i.e., for the whole
%        VineCPP toolbox) seed for random number generation.
%
% Remarks
%        * The seed is saved in the file "/private/Seed.dat".
%        * The seed is generate using the current system time.
%        * If you want to save the seed for future usage, please make a copy
%          of the file "Seed.dat". Every time random numbers are generated in
%          C++ the file will be overwritten.
%        * If the Seed should be set back to some former state, you just have
%          to replace the "Seed.dat" file in "/private/" by your copy of
%          Seed.dat you have previously generated.
%        * For further information please visit:
%           http://www.boost.org/doc/libs/1_56_0/doc/html/boost_random/reference.html
%
%
% Usage
%               SetGlobalPairCopulaParameterBounds
%
%
%
% Author: Malte Kurz

SetSeed;

end
