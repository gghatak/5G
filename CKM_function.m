function [load_BS_CKM, SINR_CKM, SNR_CKM, BS_x_CKM, BS_y_CKM, f_CKM] = CKM_function(system)

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
load_BS_CKM = zeros(1, N_BS);
%%
%Constrained K Means
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;

[index_CKM, C] = constrainedKMeans(coordinate,N_BS,N_user/(1.2*N_BS), 30);

BS_x_CKM = C(:,1);
BS_y_CKM = C(:,2);

for i = 1:N_user        
    dist_vec_CKM(i,:) = sqrt((x_vec(i) - BS_x_CKM).^2 + (y_vec(i) - BS_y_CKM).^2);
    all_BS_dist_CKM = dist_vec_CKM(i,:);
    [serving_dist_CKM(i), index_temp(i)] = min(all_BS_dist_CKM);
    interfering_dist_CKM = all_BS_dist_CKM(all_BS_dist_CKM ~= serving_dist_CKM(i));
    Interference_CKM(i) = sum(K*Ps*(sqrt(interfering_dist_CKM.^2 + H^2)).^(-alpha));
    SNR_CKM(i) = (K*Ps*sqrt(serving_dist_CKM(i)^2 + H^2)^(-alpha)) / (No); 
    SINR_CKM(i) = (K*Ps*sqrt(serving_dist_CKM(i)^2 + H^2)^(-alpha)) / (No + Interference_CKM(i));
end
for i = 1:N_user
    load_CKM(i) = 1/((1/sum(index_CKM == index_CKM(i)))*(log2(1 + SINR_CKM(i))));
    load_BS_CKM(index_CKM(i)) = load_BS_CKM(index_CKM(i)) + load_CKM(i)/sum(index_CKM == index_CKM(i));
end

[f_CKM, t_CKM] = hist(index_CKM, [1:N_BS] );
f = f_CKM/sum(f_CKM);