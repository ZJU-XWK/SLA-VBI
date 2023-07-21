function [H_eta, sigPar] = Derivation_E(sigPar, eta)
    grid = sigPar{1}.grid;
    nGrid = sigPar{1}.nGrid;
    M = sigPar{1}.M;
    N = sigPar{1}.N;
    UU = sigPar{1}.UU;
    scalar_eta = sigPar{1}.scalar_theta;
    eta_tilde = grid+scalar_eta*eta;
    % basic matrix 
    A = zeros(N, nGrid);
    for grid_index = 1:nGrid
        A(:, grid_index) = array_response(eta_tilde(grid_index), N);
    end
    F = UU*A;
    sigPar{1}.A = A;
    sigPar{1}.F = F;
    % Derivations
    grad_A_to_eta = zeros(N, nGrid);
    for grid_index = 1:nGrid
        grad_A_to_eta(:, grid_index) = 1j*pi*cos(eta_tilde(grid_index))*scalar_eta*(0:1:N-1)'.*A(:, grid_index);
    end
    grad_F_to_eta = UU*grad_A_to_eta;
    H_eta = grad_F_to_eta'*grad_F_to_eta;
    sigPar{1}.grad_F_to_theta = grad_F_to_eta;
    sigPar{1}.H_theta = H_eta;
end


