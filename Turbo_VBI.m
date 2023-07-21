function [x_VBI, sigPar] = Turbo_VBI(sigPar, Max_iter)
    y = sigPar{1}.y;
    F = sigPar{1}.F;
    [M, N] = size(F);
    PI = sigPar{1}.PI;
    x_VBI = zeros(N, Max_iter);
    % gamma parameter of noise precision 
    anoise=1e-6;
    bnoise=1e-6;
    % gamma parameter of sigal precision if zero
    a0=1*ones(N,1);
    b0=1e-5*ones(N,1);
    % gamma parameter of sigma presicion if non-zero 
    a=ones(N,1);
    b=ones(N,1);
    % initial
    PI_tilde = PI;
%     a_tilde = a.*PI+(1-PI).*a0;
%     b_tilde = b.*PI+(1-PI).*b0;
    a_tilde = sigPar{1}.a_tilde;
    b_tilde = sigPar{1}.b_tilde;
    noise_precision=sigPar{1}.noise_precision;
    for iter = 1:Max_iter
        % update for q(x)
        SigmaX = (diag(a_tilde./b_tilde)+noise_precision*(F'*F))^-1;
        MuX = SigmaX*F'*y*noise_precision;
        % update for q(gamma)
        a_tilde = PI_tilde.*a+(1-PI_tilde).*a0+1;
        b_tilde = PI_tilde.*b+(1-PI_tilde).*b0+abs(MuX).^2+(diag(SigmaX));
        % update for q(s)
        S1_BS=b.^a./gamma(a).*exp((a-1).*(psi(a_tilde)-log(b_tilde))-b.*a_tilde./b_tilde);
        S0_BS=b0.^a0./gamma(a0).*exp((a0-1).*(psi(a_tilde)-log(b_tilde))-b0.*a_tilde./b_tilde);
        S1=PI.*S1_BS;
        S2=(1-PI).*S0_BS;
        PI_tilde = S1./(S1+S2);
        % update for noise-precision
        anoise_tilde=M+anoise;
        bnoise_tilde=bnoise+norm(y-F*MuX)^2+real(trace(F*SigmaX*F'));
        noise_precision=real(anoise_tilde/bnoise_tilde);
        x_VBI(:,iter) = MuX;
    end
    sigPar{1}.x_post = MuX;
    sigPar{1}.V_post = SigmaX;
    sigPar{1}.PI = PI_tilde;
    sigPar{1}.a_tilde = a_tilde;
    sigPar{1}.b_tilde = b_tilde;
    sigPar{1}.noise_precision = noise_precision;
end

