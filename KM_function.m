function [load_BS_KM, SINR_KM, SNR_KM, BS_x_KM, BS_y_KM, f_KM] = KM_function(system)

%System Parameters
Ps = system.Ps;
Pm = system.Pm;
B = system.B;
No = system.No;
fc = system.fc;
K = system.K;
alpha = system.alpha;
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

load_BS_KM = zeros(1, N_BS);
%%
%K Means
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;
[index_KM,C] = kmeans(coordinate,N_BS);
BS_x_KM = C(:,1);
BS_y_KM = C(:,2);
for i = 1:N_user        
    dist_vec_KM(i,:) = sqrt((x_vec(i) - BS_x_KM).^2 + (y_vec(i) - BS_y_KM).^2);
    all_BS_dist_KM = dist_vec_KM(i,:);
    [serving_dist_KM(i), index_temp(i)] = min(all_BS_dist_KM);    
    interfering_dist_KM = all_BS_dist_KM(all_BS_dist_KM ~= serving_dist_KM(i));
    Interference_KM(i) = sum(K*Ps*(sqrt(interfering_dist_KM.^2 + H^2)).^(-alpha));
    SNR_KM(i) = (K*Ps*sqrt(serving_dist_KM(i)^2 + H^2)^(-alpha)) / (No); 
    SINR_KM(i) = (K*Ps*sqrt(serving_dist_KM(i)^2 + H^2)^(-alpha)) / (No + Interference_KM(i)); 
end

for i = 1:N_user
    load_KM(i) = 1/((1/sum(index_KM == index_KM(i)))*(log2(1 + SINR_KM(i))));
    load_BS_KM(index_KM(i)) = load_BS_KM(index_KM(i)) + load_KM(i)/sum(index_KM == index_KM(i));
end
[f_KM, t_KM] = hist(index_KM, [1:N_BS] );
f = f_KM/sum(f_KM);