function [xhat,xvar] = CBG_Est(rhat, rvar,message_out)
    theta = 0; %Mean of the Gaussian distribution
    var = 1; %Variance of the Gaussian distribution
    alpha = -abs(rhat-theta).^2./(var+rvar)+abs(rhat).^2./rvar;
    beta = message_out./(1-message_out).*(rvar./(var + rvar)).*exp(alpha);
    pii = 1./(1+1./beta);
    gamma = (rhat./rvar + theta./var)./(1./rvar + 1./var);
    nu = 1./(1./rvar + 1./var);
    xhat = pii.*gamma;
    xvar = pii.*(nu + abs(gamma).^2) - pii.^2 .* abs(gamma).^2;
    %xvar = pii.*nu-pii.*abs(xhat).^2.*(1./(pi*nu).*exp(-abs(gamma).^2./nu))+(1-pii).*abs(xhat).^2;
end
