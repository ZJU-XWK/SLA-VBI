function [NMSE, sigPar] = Alter_IFVBI(sigPar, Max_iter)
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
    a=ones(nGrid,1);
    b=ones(nGrid,1);
    % initial
    PI_tilde = PI;
    a_tilde = a.*PI+(1-PI).*a0;
    b_tilde = b.*PI+(1-PI).*b0;
    noise_precision = anoise/bnoise;
%     noise_precision = 1./sigPar{1}.sigma2;
    MuX=sigPar{1}.x_post;
    MuTheta = zeros(nGrid, 1);
    % Alter_VBI
    for iter = 1:Max_iter
        Ztheta = MuTheta;
        % Derivations in E-step 
        [H_theta,sigPar] = Derivation_E(sigPar, MuTheta);
        F = sigPar{1}.F;
        if (iter == 1)
            iterations = 10;
            HessX = F'*F;
            eig_x = max(eig(HessX));
        elseif (iter==2)
            iterations = 1;
            HessX = F'*F+diag(diag(H_theta).*SigmaTheta);
%             eig_x = max(eig(HessX));
        else
            iterations = 1;
            HessX = F'*F+diag(diag(H_theta).*SigmaTheta);
        end    
         % update for q(x)
        for i = 1:iterations 
            FHY = F'*y;
            Lx = zeros(20,1);
            for iter_x = 1:20
                 Lx(iter_x) = eig_x*1.01^(iter_x-1);
                 step_size = 1./(noise_precision*Lx(iter_x)+a_tilde./b_tilde);
                 MuX = step_size*noise_precision.*(FHY-HessX*MuX+Lx(iter_x).*MuX);
            end
            SigmaX = diag(1./(noise_precision*diag(HessX)+a_tilde./b_tilde));
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
            bnoise_tilde=bnoise+y'*y-2*real(MuX'*F'*y)+MuX'*HessX*MuX+trace(HessX*SigmaX);   
            noise_precision=real(anoise_tilde/bnoise_tilde);
        end
        % NMSE performance
        h_est = sigPar{1}.A*MuX;
        NMSE(iter) = norm(sigPar{1}.h - h_est, 'fro')^2/norm(sigPar{1}.h, 'fro')^2;
        % Derivations in M-step
        [Grad_theta, Hess_theta] = Derivation_M(sigPar, MuX, SigmaX); 
        % update for q(theta)
        if (iter==1)
             eig_theta = max(eig(Hess_theta));
        end
        b2 = Grad_theta+Hess_theta*Ztheta;
        Ltheta = zeros(20,1);
        for iter_theta = 1:20
        Ltheta(iter_theta) = eig_theta*1.01^(iter_theta-1);
        step_size = 1/(noise_precision*Ltheta(iter_theta)+0.5);
        MuTheta =  step_size*(noise_precision*(b2-Hess_theta*MuTheta+Ltheta(iter_theta).*MuTheta));
        end
        SigmaTheta = 1./diag(noise_precision*Hess_theta+0.5*eye(nGrid));
    end
    sigPar{1}.x_post = MuX;
    sigPar{1}.V_post = SigmaX;
    sigPar{1}.s_hat = PI_tilde;
    sigPar{1}.grid =  sigPar{1}.grid+sigPar{1}.scalar_theta*MuTheta;
end

