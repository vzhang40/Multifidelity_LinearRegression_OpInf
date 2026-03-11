function Ohat_mf = mfopinf(Xhr, Xhrdot, Xlr, Xlrdot, m1, m2, alpha, XXT)
    p = size(XXT, 1);
    r = length(alpha);
    Ohat_mf = zeros(r, p);

    for i = 1:r
    % MFMC: Hifi Data 
    indsMF = randperm(size(Xhr, 2), m2(i));
    Xh = Xhr(:, indsMF(1:m1(i)));
    Yh = Xhrdot(:, indsMF(1:m1(i)))';
    Xl1 = Xlr(:, indsMF(1:m1(i)));
    Yl1 = Xlrdot(:, indsMF(1:m1(i)))';
    Xl2 = Xlr(:, indsMF);
    Yl2 = Xlrdot(:, indsMF)';
 
    XYlm1 = (1./size(Yl1, 1))*Xl1*Yl1(:, i);
    XYlm2 = (1./size(Yl2, 1))*Xl2*Yl2(:, i);
    XYhm1 = (1./size(Yh, 1))*Xh*Yh(:, i);
    XYmulti = XYhm1 + alpha(i)*(XYlm2 - XYlm1);
    beta_multi = (XXT+1e-3*eye(p))\(XYmulti);
    Ohat_mf(i, :) = beta_multi';
    end
end

