function [Pr, C, N] = CDF_TransOT(s, t, Q, Os, pi_)
%   CDF_TransOT compute the probabilty transient operation time Q(t) 
% of a system fomulated as an CTMC Q in [0,t]. It can be presented a
% Pr(O(t)<=s)
%   The algorithm is developed [1]. The algorithm use uniformization
% technique.
% Author: Kai-Wen Tien
% Date: 06.22.2019
% Institue: Penn State University
%   s:      target value of the opertation time
%   t:      the studied transient period
%   Q:      the MxM infinitesimal generator of the system
%   Os:     the 1xM {0,1} array: Os(s) = 1 if s is operation state
%   pi_:    the 1xM initial probability of the CTMC
%
%   Pr:     the CDF Pr(O(t)<=s)
%   C:      max # of columns in the resursion for omegas
%   N:      max # of rows of in the resursion for omegas
%   revision: speed up: evaluate binomial and poisson distribution with
%   recursive function.

    % #1: normalization (randomization)
    M = size(Q,1);
    u_lam = max(-diag(Q))*1.1;
    P = Q/u_lam + eye(M);

    % #2: calculate the bound of n with a given epslion (case study setting in paper)
    epslion = 10^-4;
    N = poissinv(1-epslion/2, u_lam*t);

    % #3 calculate the bound of c and all omega values
    [C, omega, ~, ~] = recur_omega(epslion, N, P, pi_, Os);

    % get all binomal and poisson values
    binoCDF = recur_binocdf(N, C, s, t);
    poissPDF = recur_poisspdf(N, C, u_lam, t);

    % #4 sum up all probabilities
    Pr_matrix = poissPDF(:,2:end).*omega(:,2:end).*binoCDF(:,2:end);
    Pr = sum(Pr_matrix(:));
end
