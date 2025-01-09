function  [labels,centroids] = KCenters(X, K)

x_vec = X(:,1);
y_vec = X(:,2);
Xmin = min(x_vec);
Xmax = max(x_vec);
Ymin = min(y_vec);
Ymax = max(y_vec);
N_user = length(x_vec);
N_BS = K;
coordinate(:,1) = x_vec;
coordinate(:,2) = y_vec;
init_ind = randi([1, length(x_vec)]);
XBS2 = x_vec(init_ind);
YBS2 = y_vec(init_ind);
for i = 2:K    
    for i = 1:N_user        
        dist_vec_KC(i,:) = sqrt((x_vec(i) - XBS2).^2 + (y_vec(i) - YBS2).^2);
        all_BS_dist_KC = dist_vec_KC(i,:);
        [serving_dist_KC(i), index_temp(i)] = min(all_BS_dist_KC);
    end
    z = find(serving_dist_KC == max(serving_dist_KC), 1, 'first');
    XBS2 = [XBS2 x_vec(z)];
    YBS2 = [YBS2 y_vec(z)];
    clear dist_vec_KC;
end
    for i = 1:N_user        
        dist_vec_KC(i,:) = sqrt((x_vec(i) - XBS2).^2 + (y_vec(i) - YBS2).^2);
        all_BS_dist_KC = dist_vec_KC(i,:);
        [serving_dist_KC(i), labels(i)] = min(all_BS_dist_KC);
    end

centroids(:,1) = XBS2;
centroids(:,2) = YBS2;