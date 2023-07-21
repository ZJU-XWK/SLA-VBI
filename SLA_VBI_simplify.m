function [NMSE, sigPar] = SLA_VBI_simplify(sigPar, Max_iter)
    nGrid = sigPar{1}.nGrid;
    F = sigPar{1}.F;
    [M, N] = size(F);
    y = sigPar{1}.y;
    PI = sigPar{1}.PI;
    NMSE = zeros(Max_iter, 1);
    % gamma parameter of noise precision 
    anoise=1e-6;
    bnoise=1e-6;
    % gamma parameter of sigal precision if zero
    a0=1*ones(nGrid,1);
    b0=1e-5*ones(nGrid,1);
    % gamma parameter of sigma presicion if non-zero 
    a=1*ones(nGrid,1);
    b=1*ones(nGrid,1);
    % gamma parameter of off-grid parameters:theta and tao
    a_theta = 1e-6*ones(nGrid, 1);
    b_theta = 1e-6*ones(nGrid, 1);
    % initial
    PI_tilde = PI;
    a_tilde = a.*PI+(1-PI).*a0;
    b_tilde = b.*PI+(1-PI).*b0;
    noise_precision=anoise/bnoise;
%     noise_precision=1/sigPar{1}.sigma2;
    MuX=sigPar{1}.x_post;
    MuTheta = zeros(nGrid, 1);
    SigmaTheta = diag(b_theta./a_theta);
    Index_set = [];
    % Alter_VBI
    for iter = 1:Max_iter
        % Derivations in E-step
        if (iter == 1)
            [~,sigPar] = Derivation_E_simplify(sigPar, MuTheta, Index_set,true);
            HessX = F'*F;
            iterations = 10;
        else
            [H_theta,sigPar] = Derivation_E_simplify(sigPar, MuTheta, Index_set,false);
            F = sigPar{1}.F;
            HessX = F'*F;
            HessX(Index_set, Index_set) = HessX(Index_set, Index_set)+H_theta(Index_set,Index_set).*SigmaTheta(Index_set,Index_set);
            iterations = 1;
        end
         % update for q(x)
        for i = 1:iterations
            if (iter==1)
                SigmaX = (noise_precision*(HessX)+diag(a_tilde./b_tilde))^-1;
                MuX = SigmaX*noise_precision*(F'*y);
            else 
                SigmaX = (noise_precision*HessX+diag(a_tilde./b_tilde))^-1;
                MuX = SigmaX*noise_precision*(F'*y);
                %SigmaX = (noise_precision*(F'*F+H_theta.*SigmaTheta+H_tao.*SigmaTao+diag(diag(SigmaTheta))*(G_theta*F+F'*G_theta')/2+diag(diag(SigmaTao))*(G_tao*F+F'*G_tao')/2)+diag(a_tilde./b_tilde))^-1;
%                 MuX = SigmaX*noise_precision*(F'*y+0.5*diag(diag(SigmaTheta))*G_theta*y+0.5*diag(diag(SigmaTao))*G_tao*y);
            end
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
            if (iter==1)
                anoise_tilde=M+anoise;
                bnoise_tilde=bnoise+norm(y-F*MuX)^2+real(trace(F*SigmaX*F'));
            else
                anoise_tilde=M+anoise;
                bnoise_tilde=bnoise+y'*y-2*real(MuX'*F'*y)+MuX'*HessX*MuX+trace(HessX*SigmaX);   
            end
            noise_precision=real(anoise_tilde/bnoise_tilde);
        end    
        % NMSE performance
        h_est = sigPar{1}.A*MuX;
        NMSE(iter) = norm(sigPar{1}.h - h_est, 'fro')^2/norm(sigPar{1}.h, 'fro')^2;
        % Derivations in M-step
        [Grad_theta, Hess_theta] = Derivation_M(sigPar, MuX, SigmaX); 
        % update for q(theta)
        SigmaTheta_pre = (noise_precision*Hess_theta+0.5*eye(nGrid))^-1;
        MuTheta_pre =  SigmaTheta_pre*(noise_precision*(Grad_theta+Hess_theta*MuTheta));
        Index_set = find(abs(MuTheta_pre)>1e-3);
        SigmaTheta = zeros(nGrid,nGrid);
        MuTheta = zeros(nGrid,1);
        SigmaTheta(Index_set,Index_set) = SigmaTheta_pre(Index_set,Index_set);
        MuTheta(Index_set) = MuTheta_pre(Index_set);
    end
    sigPar{1}.x_post = MuX;
    sigPar{1}.V_post = SigmaX;
    sigPar{1}.s_hat = PI_tilde;
    sigPar{1}.grid =  sigPar{1}.grid+sigPar{1}.scalar_theta*MuTheta;
end

