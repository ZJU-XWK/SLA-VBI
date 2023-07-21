clc;
clear all;
% p = parpool(40);

%% System Parameters
nTxantenna = 128;
nPilot = 64;
nPath = 2;
grid_active_index = [30, 80].';
nGrid = nTxantenna;

%% Load Channel Data
load dataset/AoA_grid
load dataset/h_support
load dataset/channel_coeff_set
load dataset/AoA_set
load dataset/h_set
load dataset/UU_random_pilot
load dataset/Noise_set

%% Initialization
M = nPilot;
N = nTxantenna;
SNR_set = 0:5:20;
nSample = 1;
Ite_Max = 10;
outside_loop = 50;
Nsim = 1;
NMSE_H = zeros(outside_loop, Nsim, nSample, length(SNR_set));
AoA_GRID = zeros(nGrid,Nsim, nSample, length(SNR_set));
X_POST = zeros(nGrid, Nsim, nSample, length(SNR_set));
% AoA basis
A = zeros(N, nGrid);
for grid_index = 1:nGrid
    A(:, grid_index) = array_response(AoA_grid(grid_index, 1), N);
end
% % pilot
% permutation_S = randperm(N);
% permutation_P = randperm(N);
% S = zeros(M, N);
% P = zeros(N, N);
% p = randperm(N);
% for row_index = 1:M
%     S(row_index, permutation_S(row_index))=1;
% end    
% for row_index = 1:N
%     P(row_index,permutation_P(row_index)) = 1;
% end
% UU = S*A'*P*A';
UU = UU(1:M,:);
% noise
Noise_set = Noise_set(1:M,:);
% sensing matrix
F = UU*A;

%%
for snr_index = 3%:length(SNR_set)
    snr = SNR_set(snr_index);
    for sample_index = 1:nSample
        % channel sample
        h = h_set(:, sample_index);
        % signal
        signal = UU*h;   
        % True grid, AoA, delay
        AoA_grid_true = AoA_grid;
        AoA_grid_true(grid_active_index) = AoA_set(:, sample_index);
        for noise_index = 1:Nsim
            % noise variance
            sigma2 = (norm(signal)^2)/(M*10^(snr/10));
            % noise
            noise = sqrt(sigma2)*Noise_set(:,noise_index);
            % received signal
            y = signal+noise;
            
            %% Channel estimation
            sigPar = cell(1,1);
            sigPar{1}.M = M;
            sigPar{1}.N = N;
            sigPar{1}.UU = UU;
            sigPar{1}.PI = 0.05*ones(nGrid, 1);
            sigPar{1}.y = y;
            sigPar{1}.h = h;
            sigPar{1}.nPath = nPath;
            sigPar{1}.nGrid = nGrid;
            sigPar{1}.sigma2 = sigma2;
            sigPar{1}.noise_precision = 1;
            sigPar{1}.a_tilde = 1*ones(nGrid,1);
            sigPar{1}.b_tilde = sigPar{1}.PI+(1-sigPar{1}.PI).*1e-5.*ones(nGrid,1);
            sigPar{1}.grid = AoA_grid;
            sigPar{1}.AoA_true = AoA_set(:, sample_index);
            sigPar{1}.x_post = F'*y;
            sigPar{1}.theta_post = zeros(nGrid, 1);
            sigPar{1}.scalar_theta = pi/(N*2);
            sigPar{1}.grid_active_index = grid_active_index;
            sigPar{1}.A = A;
            sigPar{1}.F = F;
%             h = zeros(nGrid, 1);          
%             h(grid_active_index) = channel_coeff_set(:, sample_index);
%             signal1 = F*h;
%             channel esitmaiton
            [NMSE_iter, sigPar] = SLA_VBI_simplify(sigPar, outside_loop);
            NMSE_H(:,noise_index,sample_index,snr_index) = NMSE_iter;
%             for loop_index = 1:outside_loop
%                 % E-step
%                 [x_est, sigPar] = SBL(sigPar, Ite_Max);
%                 h_est = sigPar{1}.A*x_est(:,end);
%                 NMSE_iter = norm(h - h_est, 'fro')^2/norm(h, 'fro')^2; 
%                 NMSE_H(loop_index,noise_index,sample_index,snr_index) = NMSE_iter;                 
%                 % M-step
%                 [sigPar] = M_step(sigPar);
%             end
            AoA_GRID(:, noise_index,sample_index,snr_index) = sigPar{1}.grid;
            X_POST(:, noise_index,sample_index,snr_index) = sigPar{1}.x_post;
        end    
    end
end
NMSE = squeeze(mean(NMSE_H(:,:,:,:), [2, 3]));
File_name = ['result/random pilot/pilot=',num2str(nPilot),'/OGSBI/'];
% if ~isfolder(File_name)
%     mkdir(File_name);
% end
% save([File_name, 'NMSE_H'], 'NMSE_H');
% save([File_name, 'AoA_GRID'], 'AoA_GRID');
% save([File_name, 'X_POST'], 'X_POST');
%%str
% delete(p);