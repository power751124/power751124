function [binoCDF] = recur_binocdf(N, C, s, t)
    binoCDF = zeros(N+1,C+1);
    p = s/t;
    for c = 1:C
        for n = c-1:N
            if n == c-1
                binoCDF(n+1, c+1) = 1;
            else
                binoCDF(n+1, c+1) = (1-p)*binoCDF(n, c) + p* binoCDF(n,c+1);
            end
        end
    end
end