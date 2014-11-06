function [] = CheckData(u)
%CHECKDATA Checking for copula data
% Purpose:
%         The function is internally used for checking whether data input
%         is from the (d-dimensional) unit cube. It throws an error if
%         there is any real number, which is not an element of [0,1]. It
%         throws back a warning, if there is any number, which is exactly 0
%         or 1. The warning has the id TOSVC:CopulaDataOnBounds and can be
%         suppressed by the command
%         warning('off','TOSVC:CopulaDataOnBounds').
%
%
% Usage:
%           CheckData(u)
%
%
% Inputs:
%         u         = A (n x d) dimensional vector of observations /
%                     numbers, for which it should be tested whether they
%                     are copula-data (i.e., if they lie in the unit cube).
%
%
%
% Author: Malte Kurz

if sum(sum(u < 0)) || sum(sum(u > 1))
    throwAsCaller(MException('','The data has to be striclty in the intervall [0,1]'))
elseif sum(sum(u == 0)) || sum(sum(u == 1))
    warning('VineCPP:CopulaDataOnBounds','Some of the data is on the bound of the unit cube. Not all methods are stable or some functions give back non-finite results, if some input (Copula-)data is exactly 0 or 1.')
end

end
