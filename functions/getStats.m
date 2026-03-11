function bigG = getStats(Xhr, Xhrdot, Xlr, Xlrdot)
    r = size(Xhr, 1);
    bigG = cell(2, 2, r);

    for i = 1:r
    g1x = Xhr.*Xhrdot(i, :);
    g2x = Xlr.*Xlrdot(i, :);

    g1 = g1x - mean(g1x, 2);
    g2 = g2x - mean(g2x, 2);

    bigG{1, 1, i} = (g1 * g1') / (size(g1,2) - 1); % Gamma_{1, 1}
    bigG{1, 2, i} = (g1 * g2') / (size(g1,2) - 1); % Gamma_{1, 2}
    bigG{2, 1, i} = (g2 * g1') / (size(g2,2) - 1); % Gamma_{2, 1}
    bigG{2, 2, i} = (g2 * g2') / (size(g2,2) - 1); % Gamma_{2, 2}
    end
end

