function [sigPar] = M_step(sigPar)
    y = sigPar{1}.y;
    F = sigPar{1}.F;
    x_post = sigPar{1}.x_post;
    V_post = sigPar{1}.V_post;
    sigma2 = sigPar{1}.sigma2;
    grid = sigPar{1}.grid;
    nGrid = sigPar{1}.nGrid;
    A = sigPar{1}.A;
    M = sigPar{1}.M;
    N = sigPar{1}.N;
    UU = sigPar{1}.UU;
    
    %% gradient w.r.t off-grid parameters
    [~,sort_index] = sort(x_post.*conj(x_post), 'descend');
    index_amp = sort_index(1:sigPar{1}.nPath*3);
    grid_update = grid;
    grad_A_to_theta = zeros(N, nGrid);
    grad_F_to_theta = zeros(M, nGrid);
    for q = 1:length(index_amp)
        grid_index = index_amp(q);
        grad_A_to_theta(:, grid_index) = array_response(grid(grid_index, 1), N).*(0:1:N-1).'*(1j*pi*cos(grid(grid_index, 1)));
        grad_F_to_theta(:, grid_index) = UU*grad_A_to_theta(:, grid_index);
        grad_Q_to_theta = -2/sigma2*real(-x_post(grid_index)'*grad_F_to_theta(:,grid_index)'*(y-F*x_post)+grad_F_to_theta(:,grid_index)'*(F*V_post(:, grid_index)));
        grid_update(grid_index) = grid(grid_index, 1) + pi/N/40*sign(grad_Q_to_theta);
    end
    
%     % Newton update w.r.t off-grid parameters
%     [~,sort_index] = sort(x_post.*conj(x_post), 'descend');
%     index_amp = sort_index(1:sigPar{1}.nPath*2);
%     grid_update = grid;
%     grad_A_to_theta = zeros(N, nGrid);
%     grad_F_to_theta = zeros(M, nGrid);
%     hess_A_to_theta = zeros(N, nGrid);
%     hess_F_to_theta = zeros(M, nGrid);
%     grad_Q_to_theta = zeros(length(index_amp), 1);
%     hess_Q_to_theta = zeros(length(index_amp), length(index_amp));
%     for q = 1:length(index_amp)
%         grid_index = index_amp(q);
%         grad_A_to_theta(:, grid_index) = array_response(grid(grid_index, 1), N).*(0:1:N-1).'*(1j*pi*cos(grid(grid_index, 1)));
%         hess_A_to_theta(:, grid_index) = array_response(grid(grid_index, 1), N).*((1j*pi*cos(grid(grid_index, 1))*(0:1:N-1).').^2+(-1j*pi*sin(grid(grid_index, 1))*(0:1:N-1).'));
%         grad_F_to_theta(:, grid_index) = UU*grad_A_to_theta(:, grid_index);
%         hess_F_to_theta(:, grid_index) = UU*hess_A_to_theta(:, grid_index);
%         grad_Q_to_theta(q) = -2/sigma2*real(-x_post(grid_index)'*grad_F_to_theta(:,grid_index)'*(y-F*x_post)+grad_F_to_theta(:,grid_index)'*(F*V_post(:, grid_index)));
%     end
%     for q = 1:length(index_amp)
%         grid_index = index_amp(q);
%         for p = 1:length(index_amp)
%             grid_index_p = index_amp(p);
%             if (p==q)
%                 hess_Q_to_theta(q, q) = -2/sigma2*real(-x_post(grid_index)'*hess_F_to_theta(:,grid_index)'*(y-F*x_post)+x_post(grid_index)'*grad_F_to_theta(:,grid_index)'*grad_F_to_theta(:,grid_index)*x_post(grid_index))+...
%                     -2/sigma2*real(hess_F_to_theta(:,grid_index)'*(F*V_post(:, grid_index))+V_post(grid_index, grid_index)*grad_F_to_theta(:,grid_index)'*grad_F_to_theta(:,grid_index));
% %              else
% %                  hess_Q_to_theta(q, p) = -2/sigma2*real(x_post(grid_index)'*grad_F_to_theta(:,grid_index)'*grad_F_to_theta(:,grid_index_p)*x_post(grid_index_p))+...
% %                    -2/sigma2*real(V_post(grid_index_p, grid_index)*grad_F_to_theta(:,grid_index)'*grad_F_to_theta(:,grid_index_p));
%            end    
%         end    
%     end
%     grid_update(index_amp, 1) = grid(index_amp, 1) - 1./diag(hess_Q_to_theta).*grad_Q_to_theta;
    %% update sensing matrix
    A_update = A;
    for q = 1:length(index_amp)
        grid_index = index_amp(q);
        A_update(:, grid_index) = array_response(grid_update(grid_index, 1), N);
    end
    F_update = UU*A_update;
    sigPar{1}.A = A_update;
    sigPar{1}.F = F_update;
    sigPar{1}.grid = grid_update;
end

