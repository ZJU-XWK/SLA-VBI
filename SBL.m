function [x_SBL, sigPar] = SBL(sigPar,Max_iter)
    A = sigPar{1}.F;
    y = sigPar{1}.y;
    [M,N]=size(A);
    alpha=ones(N,1);
    a = 1e-6;
    b = 1e-6;
    c = 1e-6;
    d = 1e-6;
    delta = 1/sigPar{1}.sigma2;
    x_SBL = zeros(N,Max_iter);
    for iter = 1 : Max_iter
        sigma = (delta*(A'*A)+diag(alpha))^-1;
        mu = delta*sigma*A'*y;
        alpha = (1+a)./((abs(mu).^2+diag(sigma))+b);
%         delta = (c+M)/(d+((y-A*mu)'*(y-A*mu)+real(trace(A*sigma*A'))));
        x_SBL(:, iter) = mu;
    end
    sigPar{1}.x_post = mu;
    sigPar{1}.V_post = sigma;
end

