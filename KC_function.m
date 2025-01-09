function [load_BS_KC, SINR_KC, SNR_KC, BS_x_KC, BS_y_KC, f_KC] = KC_function(system)

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

load_BS_KC = zeros(1, N_BS);

%%
%Constrained K Means
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;
[index_KC, C] = KCenters(coordinate,N_BS);
BS_x_KC = C(:,1);
BS_y_KC = C(:,2);

for i = 1:N_user        
    dist_vec_KC(i,:) = sqrt((x_vec(i) - BS_x_KC).^2 + (y_vec(i) - BS_y_KC).^2);
    all_BS_dist_KC = dist_vec_KC(i,:);
    [serving_dist_KC(i), index_temp(i)] = min(all_BS_dist_KC);    
    interfering_dist_KC = all_BS_dist_KC(all_BS_dist_KC ~= serving_dist_KC(i));
    Interference_KC(i) =  sum(K*Ps*(sqrt(interfering_dist_KC.^2 + H^2)).^(-alpha));
    SNR_KC(i) = (K*Ps*sqrt(serving_dist_KC(i)^2 + H^2)^(-alpha)) / (No); 
    SINR_KC(i) = (K*Ps*sqrt(serving_dist_KC(i)^2 + H^2)^(-alpha)) / (No + Interference_KC(i)); 
end

for i = 1:N_user
    load_KC(i) = 1/((1/sum(index_KC == index_KC(i)))*(log2(1 + SINR_KC(i))));
    load_BS_KC(index_KC(i)) = load_BS_KC(index_KC(i)) + load_KC(i)/sum(index_KC == index_KC(i));
end


[f_KC, t_KC] = hist(index_KC, [1:N_BS] );
f = f_KC/sum(f_KC);