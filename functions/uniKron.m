function out = uniKron(x)
% This function takes the unique Kronecker product of a length n vector
% with itself
% 
% INPUTS: 
%   x - (n x 1) vector
% 
% OUTPUTS: 
%   out - ( n(n+1)/2 x 1 ) vector
    cnt = 1;
    n = length(x);
    xx = zeros(n*(n+1)/2,1);
    for i = 1:n
        for j = i:n
            xx(cnt) = x(i) * x(j);
            cnt = cnt + 1;
        end
    end
    out = xx;
end
