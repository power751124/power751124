function [poissPDF] = recur_poisspdf(N, C, u_lam, t)
    col1 = zeros(N+1,1);
    for n = 0:N
        % it might be too small to use exp(-u_lam*t) as the starting point
        recur_spt = 0;      % recurrent starting point
        if recur_spt <= 10^-8
            col1(n+1,1) = poisspdf(n,u_lam*t);
            recur_spt = col1(n+1,1);
        else
            col1(n+1,1) = u_lam*t/(n)*col1(n,1);
        end
    end
    cols = col1.*tril(ones(N+1,C));
    poissPDF = [col1, cols];
    
end
