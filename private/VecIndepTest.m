function [pValue,TestStat,S] = VecIndepTest(Xdata,Ydata,N)
%VECINDEPTEST Testing vectorial independence
% Purpose
%        The function performs a test on vectorial independence between
%        two vectors of random variables using the test of Quessy (2010), 
%        which is a Cramer-von Mises test statistic. As proposed by Quessy
%        (2010) approximate p-values are obtained by using the multiplier
%        bootstrap approach.
% 
%
% Usage
%             [pValue,TestStat,S] = VecIndepTest(Xdata,Ydata,N)
%
% 
% Inputs
%       Xdata                     = A (n x d1) dimensional matrix of 
%                                   observations from the first vector.
%       Ydata                     = A (n x d2) dimensional matrix of 
%                                   observations from the second vector.
%       N                         = The number of boostrapped samples,
%                                   which should be used to approximate the
%                                   p-value of the test by using the
%                                   multiplier bootstrap method.
%
% 
% Outputs
%      pValue                     = The p-value of the test.
%      TestStat                   = The value of the test statistic.
%      S                          = The boostrapped values of the test
%                                   statistics.
%
%
% References
% [1]  Kojadinovic, I. and M. Holmes (2009), "Tests of independence
%      among continuous random vectors based on Cramér-von Mises
%      functionals of the empirical copula process", Journal of
%      Multivariate Analysis 100(6), pp. 1137-1154.
% [2]  Kosorok, M. R. (2008), Introduction to Empirical Processes and
%      Semiparametric Inference, Springer Series in Statistics, New
%      York, NY: Springer.
% [3]  Kurz, M. (2013), "Tests on the partial copula", Unpublished    
%      Master's Thesis, Department of Statistics, Ludwig-Maximilians-
%      University Munich.
% [4]  Quessy, J.-F. (2010), "Applications and asymptotic power of
%      marginal-free tests of stochastic vectorial independence",
%      Journal of Statistical Planning and Inference 140(11), pp. 3058-
%      3075.
% [5]  Segers, J. (2012), "Asymptotics of empirical copula processes
%      under non-restrictive smoothness assumptions", Bernoulli 18(3),
%      pp. 764-782.
% [6]  van der Vaart, A. W. and J. A. Wellner (1996), Weak Convergence
%      and Empirical Processes -- With Applications to Statistics,
%      Springer Series in Statistics, New York [u.a.]: Springer.
% 
%
%
% Author: Malte Kurz

[n,d1] = size(Xdata);
[m,d2] = size(Ydata);

if not(n==m)
    error('The number of observations from both vectors has to be equal.')
end


U = zeros([n,d1]);
V = zeros([n,d2]);

% Transforming the data to ranks.
for i = 1:d1
    U(:,i) = tiedrank(Xdata(:,i))./(n);
end 

for i = 1:d2
    V(:,i) = tiedrank(Ydata(:,i))./(n);
end


[n,~] = size(U);
S = zeros(N,1);

[TestStat,V1,V2,V1dot,V2dot,V1dotdot,V2dotdot] = CvMTestStatCPP(U,V);

R = (V1dotdot - repmat(V1dot,n,1) - repmat(V1dot',1,n) + V1) .* (V2dotdot - repmat(V2dot,n,1) - repmat(V2dot',1,n) + V2);

%xi = randn([n,N]);
xi = RandNormal(n,N);

for k = 1:N
    S(k) = xi(:,k)'*R*xi(:,k)./n;
end

pValue = sum(S>TestStat(1))./N;

end
