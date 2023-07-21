clc;
clear all;

%% System Parameters
nTxantenna = 128;

%% Grid Parameters
nGrid = nTxantenna;
AoA_grid = zeros(nGrid,1);
sin_theta_sample=(-1+1/nGrid:2/nGrid:1-1/nGrid).';
for grid_index = 1:nGrid
    AoA_grid(grid_index) = asin(sin_theta_sample(grid_index));
end 

%% Channel parameters
nPath = 2;
grid_active_index = [30,80].';
h_support = zeros(nGrid,1);
h_support(grid_active_index) = 1;
%% Dataset Setting
nSample = 80;
Power_path = 1;
channel_coeff_set = zeros(nPath, nSample);
AoA_set = zeros(nPath, nSample);
h_set = zeros(nTxantenna, nSample);
Noise_set = sqrt(0.5)*(randn(nTxantenna,20)+1j*randn(nTxantenna,20));
for sample_index = 1:nSample
    for path_index = 1:nPath
        while abs(channel_coeff_set(path_index, sample_index))<0.5*Power_path || abs(channel_coeff_set(path_index, sample_index))>1.5*Power_path
            channel_coeff_set(path_index, sample_index) = sqrt(Power_path/2) * (randn + 1j*randn);
        end
        AoA_set(path_index, sample_index) = asin(sin(AoA_grid(grid_active_index(path_index),1))+...
            0.8*(2/nGrid*rand-1/nGrid));
        h_set(:, sample_index) = h_set(:, sample_index)+...
            channel_coeff_set(path_index, sample_index)*array_response(AoA_set(path_index, sample_index), nTxantenna);
    end    
end
%% Save Data
save('dataset/AoA_grid','AoA_grid');
save('dataset/h_support','h_support');
save('dataset/channel_coeff_set','channel_coeff_set');
save('dataset/AoA_set','AoA_set');
save('dataset/h_set','h_set');
save('dataset/Noise_set','Noise_set');




