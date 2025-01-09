function [load_BS_WKHM, SINR_WKHM, SNR_WKHM, BS_x_WKHM, BS_y_WKHM, f_WKHM] = WKHM_function(system)

%System Parameters

Ps = system.Ps;
Pm = system.Pm;
B = system.B;
No = system.No;
fc = system.fc;
K = system.K;
alpha = system.alpha;
R_vec = system.R_vec;
H = system.H;
R_I = system.R_I; 
N_u = system.N_u; 
gamma = system.gamma; 
X_min = system.X_min;
X_max = system.X_max; 
Y_min = system.Y_min;
Y_max = system.Y_max; 
N_user = system.N_user;
N_BS = system.N_BS;
x_vec = system.x_vec; 
y_vec = system.y_vec;

load_BS_WKHM = zeros(1, N_BS);

%%
%WKHM Function
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;
[index_WKHM, C] = WeightedKHarmonicMeans(coordinate,N_BS, 10);
BS_x_WKHM = C(:,1);
BS_y_WKHM = C(:,2);
for i = 1:N_user        
    dist_vec_WKHM(i,:) = sqrt((x_vec(i) - BS_x_WKHM).^2 + (y_vec(i) - BS_y_WKHM).^2);
    all_BS_dist_WKHM = dist_vec_WKHM(i,:);
    [serving_dist_WKHM(i), index_temp(i)] = min(all_BS_dist_WKHM);    
    interfering_dist_WKHM = all_BS_dist_WKHM(all_BS_dist_WKHM ~= serving_dist_WKHM(i));
    Interference_WKHM(i) = sum(K*Ps*(sqrt(interfering_dist_WKHM.^2 + H^2)).^(-alpha));

    SNR_WKHM(i) = (K*Ps*sqrt(serving_dist_WKHM(i)^2 + H^2)^(-alpha)) / (No);
    SINR_WKHM(i) = (K*Ps*sqrt(serving_dist_WKHM(i)^2 + H^2)^(-alpha)) / (No + Interference_WKHM(i)); 
end

for i = 1:N_user
    load_WKHM(i) = 1/((1/sum(index_WKHM == index_WKHM(i)))*(log2(1 + SINR_WKHM(i))));
    load_BS_WKHM(index_WKHM(i)) = load_BS_WKHM(index_WKHM(i)) + load_WKHM(i)/sum(index_WKHM == index_WKHM(i));
end
[f_WKHM, t_WKHM] = hist(index_WKHM, [1:N_BS] );
f = f_WKHM/sum(f_WKHM);