% methodIC5.m
% Andrew J. Wiebe, 24 Jul 2019
% After Tarpanelli et al. (2012, JofH)
%
% Code for GNU Octave (Eaton et al., 2018).
%
% Objective: Follow the steps in Tarpanelli et al (2012) to use the Iman
% and Conover (1982) method to generate random rainfall time series with a
% desired level of Pearson (and Spearman) correlation.
%
% X = m by n matrix (n random timeseries of length m)
% C* = desired Pearson correlation matrix, calculated from (external) base matrix
%
% The correlation matrix R of the columns of X is calculated within this function for convenience. 
%
% References:
% Eaton, J.W., Bateman, D., Hauberg, S., Wehbring, R., 2018. GNU Octave. Edition 5 for Octave version 5.1.0. Manual for
%    high-level interactive language for numerical computations. https://www.gnu.org/software/octave/download.html. February 2019.
% Iman, R.L., Conover, W.J., 1982. A distribution-free approach to inducing rank correlation among input variables.
%    Commun. Statist.-Simula. Computa. 11(3), 311-334. https://doi.org/10.1080/03610918208812265.
% Tarpanelli, A., Franchini, M., Brocca, L., Camici, S., Melone, F., Moramarco, T., 2012. A simple approach for stochastic generation
%    of spatial rainfall patterns. J. Hydrol 472-473, 63-76. https://doi.org/10.1016/j.jhydrol.2012.09.010.
%

function [X_star, pear, spear] = methodIC5(X, C_star);

%%% calculate pearson coefficient matrix for random input matrix X
n = length(C_star);
R = zeros(n,n);

for i = 1:n
	for j = i:n
        R(i,j) = corr(X(:,i), X(:,j));
        R(j,i) = R(i,j);
	end
end

Q = chol(R,'lower');

P = chol(C_star, 'lower');

X1 = X * ((P * inv(Q))');

%%% Now rearrange columns of X to have the same ranks as the columns of X1?
[X1_sort, idx1] = sort(X1,1, 'descend'); % sort the columns of X1

X_sort = sort(X,1,'descend');

for col = 1:n
    for row = 1:length(X)
        X_star(idx1(row,col),col) = X_sort(row,col);
    end
end

%%% check overall cross-correlation for random time series

for i = 1:n
	for j = i:n
        pear(i,j) = corr(X_star(:,i), X_star(:,j));
        pear(j,i) = pear(i,j);
        spear(i,j) = spearman(X_star(:,i), X_star(:,j));
        spear(j,i) = spear(i,j);
	end
end