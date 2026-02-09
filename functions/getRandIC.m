function param = getRandIC(p, q)
% INPUTS: 
%   p - number of random amplitude parameters
%   q - number of random frequency parameters
% OUTPUTS: 
%   param - (2 x pq) parameters [a; b] assoicated with each intial
%       condition
    if q ~= 0
    % Generating random amplitudes
    a = rand(p, 1);

    % Generating random frequencies
    b = randi(5, [q, 1]); %

    % parameters
    param = [repmat(a', 1, q); repelem(b', 1, p)];
    inds = randperm(p*q);
    param = param(:, inds);
    else 
        a = rand(p, 1);
        param = [a'; zeros(1, p)];
    end
end

