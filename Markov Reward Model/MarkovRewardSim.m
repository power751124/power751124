function cum_rs = MarkovRewardSim(t, Q, R, a, iter)
%MarkovRewardSim is a simulator to simulate the cumulative rewards of Markov Reward Model 
%   t   : the time window [0,t]
%   Q   : the transient rate matrix NxN
%   R   : the reward matrix Nx1
%   a   : the initial probability of each state
%   iter: the iterations of the simulation
%%
    % determine the probablity of transient
    SR = -diag(Q);
    P = Q./SR + eye(size(Q));
    cum_rs = [];
    %
    for ii = 1:iter
        clock = 0;
        cum_r = 0;
        % detemine the inital state
        state = find(rand<cumsum(a),1,'first');
        while clock <= t
            tis =  exprnd(1/SR(state));
            if clock + tis <= t
                cum_r = cum_r + R(state)*tis;
            else
               cum_r = cum_r + R(state)*(t - clock);
            end
            state = find(rand<cumsum(P(state,:)),1,'first');
            clock = clock + tis;
        end
        cum_rs = [cum_rs;cum_r];
    end

end

