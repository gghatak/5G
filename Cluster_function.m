function [No_cluster, Cluster, subgraphs] = Cluster_function(system, threshold_users, algo)

system.threshold_users = threshold_users;

switch algo
    case 1
        [load_algo, SINR_algo, SNR_algo, BS_x_algo, BS_y_algo, f_algo] = KM_function(system); 
    case 2
         [load_algo, SINR_algo, SNR_algo, BS_x_algo, BS_y_algo, f_algo] = CKM_function(system); 
    case 3
         [load_algo, SINR_algo, SNR_algo, BS_x_algo, BS_y_algo, f_algo] = KHM_function(system);
    case 4
         [load_algo, SINR_algo, SNR_algo, BS_x_algo, BS_y_algo, f_algo] = WKHM_function(system);    
    case 5
         [load_algo, SINR_algo, SNR_algo, BS_x_algo, BS_y_algo, f_algo] = KC_function(system); 
end
        
  for i = 1: length(BS_x_algo)
            for j = 1:length(BS_x_algo)
                Gain_Matrix_algo(i,j) = (((BS_x_algo(i) - BS_x_algo(j)))^2 + (BS_y_algo(i) - BS_y_algo(j))^2)^-2;
            end
  end
 
  Graph_nodes_algo = [1:length(BS_x_algo)];
  [subgraphs] = k_flex_clustering(Graph_nodes_algo, Gain_Matrix_algo, f_algo, system.threshold_users);    

%New SNR Calculation
  for i = 1:system.N_user
    dist_vec_algo(i,:) = sqrt((system.x_vec(i) - BS_x_algo).^2 + (system.y_vec(i) - BS_y_algo).^2);
    all_BS_dist_algo = dist_vec_algo(i,:);    
    [serving_dist_algo(i), index_temp(i)] = min(all_BS_dist_algo);
        for j = 1:size(subgraphs,1)
            for k = 1:size(subgraphs,2)
                    if subgraphs(j,k) ==index_temp(i)
                        Clustered_serving_distances = sqrt((system.x_vec(i) - BS_x_algo(subgraphs(j,subgraphs(j,:)~=0))).^2 + (system.y_vec(i) - BS_y_algo(subgraphs(j,subgraphs(j,:)~=0))).^2);
                        Clustered_power(i) = sum(system.K*system.Ps*(sqrt(Clustered_serving_distances.^2 + system.H^2).^(-system.alpha)));
                    end
            end
        end
        Total_power(i) = sum(system.K*system.Ps*(sqrt(all_BS_dist_algo.^2 + system.H^2)).^(-system.alpha));
        Clustered_interference(i) = Total_power(i) - Clustered_power(i);
        Clustered_SNR_algo(i) = Clustered_power(i) / (system.No);
        Clustered_SINR_algo(i) = Clustered_power(i) / (system.No + Clustered_interference(i)); 
  end    

        No_cluster = min(pow2db(SNR_algo*system.No*30e3/system.B)) + 30;
        Cluster = min(pow2db(Clustered_SNR_algo*system.No*30e3/system.B)) + 30;
end