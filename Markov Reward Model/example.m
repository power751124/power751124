%% example #1
clc;clear;
Q = [-1.3, 1.3, 0, 0;
     0, -1.2, 1.2, 0;
     0, 0, -.5, .5;
     .5, 0, 0, -.5];
Os = [1,0,0,0];
pi_ = [1,0,0,0];
t = 10;
s = .5*t;
%% Simulation
cum_r = MarkovRewardSim(t,Q,Os,pi_,10000);
pr = sum(cum_r<=s)/sum(ones(size(cum_r)));
cdfplot(cum_r)
lb = min(cum_r);
ub = max(cum_r);
xlabel('capacity in [0, t]')
ylabel('probability')
title('')
%% Analytical Model
s = linspace(lb,ub,100);
for ii = 1:length(s)
    Pr(ii) = CDF_TransOT(s(ii), t, Q, Os, pi_);
end
hold on
scatter(s,Pr,'x');
legend('simu','approx')
hold off
