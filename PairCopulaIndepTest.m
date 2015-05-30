function [pValue,TestStat] = PairCopulaIndepTest(u1,u2)
%PAIRCOPULAINDEPTEST Testing on the independence pair-copula
% Purpose
%        The function performes the bivariate independence test for copula
%        data of Genest and Favre (2007).
%
%
% Usage
%       [pValue,TestStat] = PairCopulaIndepTest(u1,u2)
%
%
% Inputs
%       u1                = A (n x 1) dimensional vector of values lying in [0,1].
%       u2                = A (n x 1) dimensional vector of values lying in [0,1].
%
%
% Outputs
%      pValue             = The p-value for the independence test.
%      TestStat           = The value of the test statistic.
%
%
% References
%      Genest, C. and A. Favre (2007), "Everything You Always Wanted to Know
%      about Copula Modeling but Were Afraid to Ask", Journal of Hydrologic
%      Engineering 12(4), pp. 347-368.
%
%
%
% Author: Malte Kurz

% Check the (Copula-)data input.
CheckData([u1,u2])

[pValue,TestStat] = VineCopulaMatlab(5,u1,u2);

end
