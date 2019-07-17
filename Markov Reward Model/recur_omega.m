function [C, omega, err, omega_C] = recur_omega(epslion, N, P, pi_, Os)
%   recur_omega recursion function of omega in the algorithm developed in [1]. 
% The algorithm use uniformization technique to compute transient operational time.
%
% [1] De Souza, E. Silva. "Calculating cumulative operational time distributions of 
% repairable computer systems." IEEE Transactions on Computers 100.4 (1986): 322-332.
%
% Author: Kai-Wen Tien
% Date: 06.22.2019
% Institue: Penn State University
%   epslion:    the bound of n
%   N:          max # of rows of in the resursion for omegas
%   P:          uniformization probability matrix of Q
%   pi_:        the 1xM initial probability of CTMC
%   Os:         the 1xM {0,1} array: Os(s) = 1 if s is operation state
%
%   C:          max # of columns in the resursion for omegas
%   omega:      the omega values in algorithm [1]
%   err:        the estimated error
%   omega_C:    the 3d omega matrix:(N+1) x (C+1) x M
%   revision 1.0
    
    omega_C = [];
    M = length(pi_);    % size of the markov model
    C = 0;              % initial C = 0
    omega_c = zeros(N+1, M);
    for n = 0:N
        if n == 0
            omega_c(n+1,:) = pi_.*(Os==1);
        else
            omega_o = omega_c(n,:)*P(:,Os==1);
            omega_f = zeros(1,M)*P(:,Os==0);
            omega_c(n+1,Os==1) = omega_o;
            omega_c(n+1,Os==0) = omega_f;
        end
    end
    omega_C(:,1,:) = omega_c;
    err = 1 - sum(sum(omega_C(end,:,:)));
    
    % generate C > 0
    while err >= epslion/2 && C <= N+1
        C = C+1;
        omega_c = zeros(N+1, M);
        for n = C-1:N
            if n == 0
                omega_c(n+1,:) = pi_.*(Os==0);
            else
                omega_o = omega_c(n,:)*P(:,Os==1);
                omega_f = squeeze(omega_C(n,C,:))'*P(:,Os==0);
                omega_c(n+1,Os==1) = omega_o;
                omega_c(n+1,Os==0) = omega_f;
            end
        end
        omega_C(:,C+1,:) = omega_c;
        err = 1-sum(sum(omega_C(end,:,:)));
    end
    omega = sum(omega_C,3);
end
