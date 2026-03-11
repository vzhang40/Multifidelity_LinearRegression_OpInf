function [m1, m2, a1, a2] = mflrOptimalSampleSize(bigG, w, p)
% This is a bifidelity example
r = size(bigG, 3);
m1 = zeros(r, length(p));
m2 = zeros(r, length(p));
a1 = zeros(r, 1);
a2 = zeros(r, 1);
for i = 1:r
    A2 = bigG{1, 2, i}*pinv(bigG{2, 2, i});
    a1(i) = trace(bigG{1, 1, i} - A2*bigG{2, 1, i});
    a2(i) = trace(A2*bigG{2, 1, i});
    denom = (sqrt(a1(i).*w(1)) + sqrt(a2(i).*w(2)));
    m1(i, :) = floor( p*(sqrt(a1(i)./w(1)) / denom) );
    m2(i, :) = floor( p*(sqrt(a2(i)./w(2)) / denom) );
end
m1(m1 < 1) = 1;
m2(m2 - m1 < 0) = m1(m2 - m1 < 0) + 1;
end

