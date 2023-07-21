function [NMSE, sigPar] = OGSBI(sigPar, Max_iter)
    nGrid = sigPar{1}.nGrid;
    F = sigPar{1}.F;
    [M, N] = size(F);
    y = sigPar{1}.y;
    % gamma parameter of noise precision
    anoise=1e-6;
    bnoise=1e-6;
    % gamma parameter of sigma presicion
    alpha = ones(nGrid, 1);
    a=1*ones(nGrid,1);
    b=1e-2*ones(nGrid,1);
    % initial
    noise_precision = anoise/bnoise;
%     noise_precision=1/sigPar{1}.sigma2;
    MuX=sigPar{1}.x_post;
    MuTheta = zeros(nGrid, 1); 
    [~, sigPar] = Derivation_E(sigPar, MuTheta);
    A = sigPar{1}.F;
    B = sigPar{1}.grad_F_to_theta;
    idx = [];
     % OGSBI
    for iter = 1:Max_iter
        % Derivations in E-step
        F(:,idx) = A(:,idx) + B(:,idx) * diag(MuTheta(idx));
         % update for q(x)
        SigmaX = (noise_precision*(F'*F)+diag(alpha))^-1;
        MuX = SigmaX*noise_precision*(F'*y);
        % update for q(gamma)
        alpha = (2*b)./(sqrt(1+4*b.*(abs(MuX).^2+diag(SigmaX)))-1);
%        alpha = (1+a)./(abs(MuX).^2+diag(SigmaX)+b);
        % update for noise-precision
        anoise_tilde=M+anoise;
%         gamma1 = 1 - real(diag(SigmaX)).*alpha;
        bnoise_tilde=bnoise+norm(y-F*MuX)^2+real(trace(F*SigmaX*F'));
%         bnoise_tilde = bnoise+norm(y-F*MuX)^2+ 1 / noise_precision * sum(gamma1);
        noise_precision=real(anoise_tilde/bnoise_tilde);
%       sigPar = OGSBI_Sen_Mat(sigPar, MuTheta);
        % Derivations in M-step
        [~, idx] = sort(abs(MuX), 'descend');
        idx = idx(1:sigPar{1}.nPath);
        P = real(conj(B(:,idx)'*B(:,idx)).*(MuX(idx)*MuX(idx)'+SigmaX(idx,idx)));
        v = real(conj(MuX(idx)) .* (B(:,idx)' * (y - A * MuX)))- real(diag(B(:,idx)'*A*SigmaX(:,idx)));
        % update for q(theta)
%         SigmaTheta = (noise_precision*Hess_theta)^-1;
%         MuTheta =  SigmaTheta*(noise_precision*(Grad_theta+Hess_theta*MuTheta));
        MuTheta = zeros(N, 1);
        MuTheta(idx) = P^-1*v;
        MuTheta = max(MuTheta, -1);
        MuTheta = min(MuTheta, 1);
    end
    % SBL estimation
    sigPar = OGSBI_Sen_Mat(sigPar, MuTheta);
    sigPar{1}.F = sigPar{1}.UU*sigPar{1}.A;
    Ite_Max = 20;
    [x_est, sigPar] = SBL(sigPar, Ite_Max);
    % NMSE performance
    NMSE_iter = zeros(Ite_Max, 1);
    for ite = 1:Ite_Max
        h_est = sigPar{1}.A*x_est(:,ite);
        NMSE_iter(ite) = norm(sigPar{1}.h - h_est, 'fro')^2/norm(sigPar{1}.h, 'fro')^2;
    end
    NMSE = min(NMSE_iter);
    % save data
    sigPar{1}.x_post = MuX;
    sigPar{1}.V_post = SigmaX;
    grid = sigPar{1}.grid;
    sigPar{1}.grid =  grid+sigPar{1}.scalar_theta*MuTheta;
end

