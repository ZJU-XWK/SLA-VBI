function [x_est,sigPar] = Turbo_CS(sigPar,Ite_Max)
    y = sigPar{1}.y;
    F = sigPar{1}.F;
    M = sigPar{1}.M;
    nGrid = sigPar{1}.nGrid;
    PI = sigPar{1}.PI;   
    %% initialize
    x_pri = zeros(nGrid, 1);
    V_pri = 2*ones(nGrid, 1);
    anoise = 1e-6;
    bnoise = 1e-6;
    noise_precision = sigPar{1}.noise_precision;
    sigma2 = 1/noise_precision;
%     sigma2 = sigPar{1}.sigma2;
%     NMSEA = zeros(Ite_Max+1, 1);
%     NMSEB = zeros(Ite_Max, 1);
    x_est = zeros(nGrid, Ite_Max+1);
    %% Turbo CS Algortihm
    for iter = 1:Ite_Max
       %% Module A: LMMSE
        V_post = (F'*F/sigma2+diag(1./V_pri))^-1;
        x_post = V_post*(x_pri./V_pri+F'*y/sigma2);
        % diagnal approximation
        V_post = diag(V_post); 
        % extrinsic message
        V_ext = 1./(1./V_post-1./V_pri);
        x_ext = V_ext.*(x_post./V_post-x_pri./V_pri);
        x_pri = x_ext;  
        V_pri = V_ext;
        x_est(:,iter) = x_post;
        % NMSE performance of Module A  
%         H_post = zeros(M, N);
%         for n_index = 1:N
%             H_post(:,n_index) = sigPar{1}.A*(diag(sigPar{1}.D(n_index, :))*x_post);
%         end
%         NMSEA(iter) = norm(sigPar{1}.H - H_post, 'fro')^2/norm(sigPar{1}.H, 'fro')^2; 
        %% Module B: Message Passing
        message_in_s = message_in(x_pri, V_pri);
        message_out_s = PI;

        [x_post,V_post] = CBG_Est(x_pri, V_pri,message_out_s);
        err_index = find(V_post >= V_pri);
        V_ext = (1./V_post-1./V_pri).^-1;
        x_ext = V_ext./V_post .* x_post -V_ext./V_pri .* x_pri;
        V_ext(err_index) = 2 * V_pri(err_index);
        x_ext(err_index) = 0;
        x_pri = x_ext;
        V_pri = V_ext;
%         NMSE performance of Module B  
%         H_post = zeros(M, N);
%         for n_index = 1:N
%             H_post(:,n_index) = sigPar{1}.A*(diag(sigPar{1}.D(n_index, :))*x_post);
%         end   
%         NMSEB(iter) = norm(sigPar{1}.H - H_post, 'fro')^2/norm(sigPar{1}.H, 'fro')^2;
    end
    % update for noise_precision
    anoise_tilde = anoise+M;
    bnoise_tilde = bnoise+norm(y-F*x_post)^2+real(trace(F*diag(V_post)*F'));
    noise_precision = real(anoise_tilde/bnoise_tilde);
    sigma2 = 1/noise_precision;
   %% Module A: LMMSE 
    V_post = (F'*F/sigma2+diag(1./V_pri))^-1;
    x_post = V_post*(x_pri./V_pri+F'*y/sigma2);
    x_est(:, end) = x_post;
    sigPar{1}.x_post = x_post;
    sigPar{1}.V_post = V_post;
    s_post = real(message_in_s.*message_out_s./(message_in_s.*message_out_s+(1-message_in_s).*(1-message_out_s)));    
    sigPar{1}.s_post = s_post;
    sigPar{1}.noise_precision = noise_precision;
end


