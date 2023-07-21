function [sigPar] = OGSBI_Sen_Mat(sigPar, theta)
    grid = sigPar{1}.grid;
    nGrid = sigPar{1}.nGrid;
    N = sigPar{1}.N;
    scalar_eta = sigPar{1}.scalar_theta;
    theta_tilde = grid+scalar_eta*theta;
    % basic matrix 
    A = zeros(N, nGrid);
    for grid_index = 1:nGrid
        A(:, grid_index) = array_response(theta_tilde(grid_index), N);
    end
    sigPar{1}.A = A;
end


