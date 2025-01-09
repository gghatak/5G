function [load_BS_KHM, SINR_KHM, SNR_KHM, BS_x_KHM, BS_y_KHM, f_KHM] = KHM_function(system)

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

load_BS_KHM = zeros(1, N_BS);

%%
%K Means
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;
[index_KHM, C] = KHarmonicMeans(coordinate,N_BS, 25);
BS_x_KHM = C(:,1);
BS_y_KHM = C(:,2);

for i = 1:N_user        
    dist_vec_KHM(i,:) = sqrt((x_vec(i) - BS_x_KHM).^2 + (y_vec(i) - BS_y_KHM).^2);
    all_BS_dist_KHM = dist_vec_KHM(i,:);
    [serving_dist_KHM(i), index_temp(i)] = min(all_BS_dist_KHM);  
    interfering_dist_KHM = all_BS_dist_KHM(all_BS_dist_KHM ~= serving_dist_KHM(i));
    Interference_KHM(i) = sum(K*Ps*(sqrt(interfering_dist_KHM.^2 + H^2)).^(-alpha));
    SNR_KHM(i) = (K*Ps*sqrt(serving_dist_KHM(i)^2 + H^2)^(-alpha)) / (No);
    SINR_KHM(i) = (K*Ps*sqrt(serving_dist_KHM(i)^2 + H^2)^(-alpha)) / (No + Interference_KHM(i)); 
end

for i = 1:N_user
    load_KHM(i) = 1/((1/sum(index_KHM == index_KHM(i)))*(log2(1 + SINR_KHM(i))));
    load_BS_KHM(index_KHM(i)) = load_BS_KHM(index_KHM(i)) + load_KHM(i)/sum(index_KHM == index_KHM(i));
end

[f_KHM, t_KHM] = hist(index_KHM, [1:N_BS] );
f = f_KHM/sum(f_KHM);