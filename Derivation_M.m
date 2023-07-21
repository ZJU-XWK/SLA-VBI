function [Grad_theta, Hess_theta] = Derivation_M(sigPar, MuX, SigmaX)
    y = sigPar{1}.y;
    F = sigPar{1}.F;
    grad_F_to_theta = sigPar{1}.grad_F_to_theta;
    H_theta = sigPar{1}.H_theta;
    Grad_theta = real(diag(MuX)'*(grad_F_to_theta'*(y-F*MuX))-diag(grad_F_to_theta'*F*SigmaX));
    Hess_theta = real((MuX*MuX'+SigmaX).'.*H_theta);
end