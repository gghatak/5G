clc; clear all; 

%%
% Description %

% This is the master code that calls different placement algorithms.
% Inputs: Number of radio nodes (RNs), the deployment area, and  the number of network realizations.
% Output: Locations of the RNs, minimum and mean SINR, SNR, RSRP across all users and network
% realizations and relative user association and cell load across RNs.
% The system parameters are defined in the structure "system". 
% The placement algorithm functions have the format
% <algo_name>_function(system). 5 different algo_name are used: k-Means
% (KM),k-harmonic means (KHM), k-centers (KC), weighted KHM (WKHM), and
% constrained KM (CKM). Each of these functions return a vector of the
% form: [load_BS_<algo_name>, SINR_<algo_name>, SNR_<algo_name>,
% BS_x_<algo_name>, BS_y_<algo_name>, f_<algo_name>]. They are described
% below:
% load_BS_<algo_name> represents the load per RN. This takes into account
% not only the number of users assocaited to the RNs but also their relative
% locations. The size of this output is the same as the number of RNs.
% SINR_<algo_name> and SNR_<algo_name> represent the SINR and the SNR of
% the users. The size of these outputs are the same as the number of users.
% BS_x_<algo_name> and BS_y_<algo_name> respectively represent the
% x-coodinate and the y-coordinate of the BSs.
% Finally, f_<algo_name> represents the number of users that are connected
% to each RNs. The size of this output is same as the number of RNs.

%%
%Inputs
system.N_BS = 10; % Number of radio nodes (RNs)
system.N_user = 100; % Number of users
iterations = 10; % Number of network realizations to average out
Deployment_Area = 10000;


%%
% System Parameters defined as a structure

%Approximation as a square grid
system.X_min = -sqrt(Deployment_Area)/2;
system.X_max = sqrt(Deployment_Area)/2;
system.Y_min = -sqrt(Deployment_Area)/2;
system.Y_max = sqrt(Deployment_Area)/2;

system.Ps = db2pow(23-30); %Transmit power (23dbm) in Watts for indoor nodes
system.Pm = db2pow(40-30); %Transmit power (40dbm) in Watts for outdoor Macro cells
system.B = 100e6; %Total bandwidth per cell
system.No = db2pow(-176 - 30 + 9)*system.B; %Noise density -176 dbm per Hz
system.fc = 3.2e9; %Carrier frequency

%Propagation Parameters
system.K = (3e8/(4*pi*system.fc))^2; %Path loss constant for indoor propagation (3GPP)
system.alpha = 3; %Path loss exponent for indoor propagation (3GPP) wthout wall material details

%Dimensions and Traffic Requirements
system.R_vec = linspace(5,15,30); %Radius of deployment area
system.H = 4; %Height of the floor
system.R_I = 20; %Interference range
system.N_u = 20; %Number of users
system.gamma = db2pow(-8); %SINR Threshold for successful decoding of control channel


%%
% Generation of skewed user locations
optionx = 2; %Model/distribution of skewness of the x-coordinate
optiony = 3; %Model/distribution of skewness of the y-coordinate

%The function generate_user_locations return the locations of N_user users
%skewed as per optionx and optiony
[system.x_vec, system.y_vec] =  generate_user_locations(system.X_max, system.Y_max, system.X_min, system.Y_min, system.N_user, optionx, optiony);

%%
%Main body of the code
for iter = 1:iterations

    disp('Running KM')
    [load_BS_KM, SINR_KM, SNR_KM, BS_x_KM, BS_y_KM, f_KM] = KM_function(system);

    disp('Running CKM')
    [load_BS_CKM, SINR_CKM, SNR_CKM, BS_x_CKM, BS_y_CKM, f_CKM] = CKM_function(system);

    disp('Running KHM')
    [load_BS_KHM, SINR_KHM, SNR_KHM, BS_x_KHM, BS_y_KHM, f_KHM] = KHM_function(system);

    disp('Running WKHM')
    [load_BS_WKHM, SINR_WKHM, SNR_WKHM, BS_x_WKHM, BS_y_WKHM, f_WKHM] = WKHM_function(system);

    disp('Running KC')
    [load_BS_KC, SINR_KC, SNR_KC, BS_x_KC, BS_y_KC, f_KC] = KC_function(system);

%%
    % Calculation of relevant metrics in vectorized format. Each metric contains
    % 5 elements, each cooresponding to each of the algorithms.

    %Minimum SINR across all the users in one network realization
    Min_SINR(iter,:) = [min(pow2db(SINR_KM)) min(pow2db(SINR_CKM)) min(pow2db(SINR_KHM)) min(pow2db(SINR_WKHM)) min(pow2db(SINR_KC))];
    
    %Minimum SNR across all the users in one network realization
    Min_SNR(iter,:) = [min(pow2db(SNR_KM)) min(pow2db(SNR_CKM)) min(pow2db(SNR_KHM)) min(pow2db(SNR_WKHM)) min(pow2db(SNR_KC))];
    
    %Mean SINR across all the users in one network realization
    Mean_SINR(iter,:) = [mean(pow2db(SINR_KM)) mean(pow2db(SINR_CKM)) mean(pow2db(SINR_KHM)) mean(pow2db(SINR_WKHM)) mean(pow2db(SINR_KC))];
    
    %Mean SNR across all the users in one network realization
    Mean_SNR(iter,:) = [mean(pow2db(SNR_KM)) mean(pow2db(SNR_CKM)) mean(pow2db(SNR_KHM)) mean(pow2db(SNR_WKHM)) mean(pow2db(SNR_KC))];
    
    %Mean load across all the RNs in one network realization
    Mean_Load(iter,:) = [mean(load_BS_KM) mean(load_BS_CKM) mean(load_BS_KHM) mean(load_BS_WKHM) mean(load_BS_KC)];
    
    %Maximum load across all the RNs in one network realization
    Max_Load(iter,:) = [max(load_BS_KM) max(load_BS_CKM) max(load_BS_KHM) max(load_BS_WKHM) max(load_BS_KC)];
    
    %Use minimum SNR to calculate the minimum RSRP. Assume 30 kHz has the tone
    %width
    RSRP(iter,:) = [min(SNR_KM)*system.No*30e3/system.B min(SNR_CKM)*system.No*30e3/system.B min(SNR_KHM)*system.No*30e3/system.B min(SNR_WKHM)*system.No*30e3/system.B min(SNR_KC)*system.No*30e3/system.B];
    RSRP_vec(iter, :) = pow2db(RSRP(iter,:))+ 30; %Convert to dBm
    
    
    %Use mean SNR to calculate the mean RSRP. Assume 30 kHz has the tone
    %width
    Mean_RSRP(iter,:) =  [mean(SNR_KM)*system.No*30e3/system.B mean(SNR_CKM)*system.No*30e3/system.B mean(SNR_KHM)*system.No*30e3/system.B mean(SNR_WKHM)*system.No*30e3/system.B mean(SNR_KC)*system.No*30e3/system.B];
    Mean_RSRP_vec(iter, :) = pow2db(Mean_RSRP(iter,:)) + 30; %Convert to dBm
    
    
    %The ratio of the number of users connected to the most admitted RN to the
    %least admitted RN
    Associate(iter,:) = [max(f_KM)/min(f_KM) max(f_CKM)/min(f_CKM) max(f_KHM)/min(f_KHM) max(f_WKHM)/min(f_WKHM) max(f_KC)/min(f_KC)];
end

%%
%Displaying the results

%Vectors collecting relevant metrics together
Signals = [mean(Min_SINR); mean(Min_SNR); mean(Mean_SINR); mean(Mean_SNR); mean(RSRP_vec); mean(Mean_RSRP_vec)]; %These correspond to the signal strengths and interference - depend only on RN numbers and locations
Loads = [mean(Max_Load) ; mean(Mean_Load); mean(Associate)]; %These depend on number of users and their locations

figure
bar(Signals)
legend('KM', 'CKM', 'KHM', 'WKHM', 'KC')
set(gca, 'XTickLabel', {'Min SINR' 'Min SNR' 'Mean SINR' 'Mean SNR' 'Min RSRP' 'Mean RSRP'})

figure
bar(Loads)
legend('KM', 'CKM', 'KHM', 'WKHM', 'KC')
set(gca, 'XTickLabel', {'Max Load' 'Mean SNR' 'Association Ratio'})