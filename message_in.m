function [message_in] = message_in(rhat,rvar)
    theta = 0;%Mean of the Gaussian distribution
    var = 1;%Variance of the Gaussian distribution
    tmpVar = rvar + var;
    a = sqrt(tmpVar./rvar).*exp(abs(rhat-theta).^2./tmpVar-abs(rhat).^2./rvar); 
    message_in = 1./(1+a); 
end