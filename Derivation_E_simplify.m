function [H_theta, sigPar] = Derivation_E_simplify(sigPar, theta, Index_set, init)
    grid = sigPar{1}.grid;
    nGrid = sigPar{1}.nGrid;
    N = sigPar{1}.N;
    UU = sigPar{1}.UU;
    scalar_theta = sigPar{1}.scalar_theta;
    theta_tilde = grid+scalar_theta*theta;
    % basic matrix 
    A = sigPar{1}.A;
    F = sigPar{1}.F;
    if (init)
        grad_A_to_theta = zeros(N, nGrid);
        for grid_index = 1:nGrid
            grad_A_to_theta(:, grid_index) = 1j*pi*cos(theta_tilde(grid_index))*scalar_theta*(0:1:N-1)'.*A(:, grid_index);
        end
        grad_F_to_theta = UU*grad_A_to_theta;
        H_theta = grad_F_to_theta'*grad_F_to_theta;
        sigPar{1}.grad_A_to_theta_init = grad_A_to_theta;
        sigPar{1}.grad_F_to_theta_init = grad_F_to_theta;
    else
        grad_A_to_theta = sigPar{1}.grad_A_to_theta_init;
        grad_F_to_theta = sigPar{1}.grad_F_to_theta_init;
        for index = 1:length(Index_set)
            grid_index = Index_set(index);
            A(:, grid_index) = array_response(theta_tilde(grid_index), N);
            F(:, grid_index) = UU*A(:, grid_index);
            grad_A_to_theta(:, grid_index) = 1j*pi*cos(theta_tilde(grid_index))*scalar_theta*(0:1:N-1)'.*A(:, grid_index);
        end
        grad_F_to_theta(:, Index_set) = UU*grad_A_to_theta(:, Index_set);
        H_theta = grad_F_to_theta'*grad_F_to_theta;
    end    
    
    % save 
    sigPar{1}.A = A;
    sigPar{1}.F = F;
    sigPar{1}.grad_A_to_theta = grad_A_to_theta;
    sigPar{1}.grad_F_to_theta = grad_F_to_theta;
    sigPar{1}.H_theta = H_theta;
end


