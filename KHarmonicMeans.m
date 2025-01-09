function  [labels,centroids] = KHarmonicMeans(X, K, maxiter)

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
[~,C] = kmeans(coordinate,N_BS);
XBS2 = C(:,1);
YBS2 = C(:,2);
count = 0;
p = 2;
while count < maxiter
            %Clustering Update
             for i = 1:N_user
                 dist_vec2(i,:) = max(1, sqrt((x_vec(i) - XBS2).^2 + (y_vec(i) - YBS2).^2));
                 [min_dist(i), ind] = min(dist_vec2(i,:));        
                 labels(i) = ind;
             end
            %Centroid Update
            for j = 1:N_BS 
                for i = 1:N_user
                        Numx(i) = x_vec(i) * (1 / ( (dist_vec2(i,j))^(p+2)  * sum(1./(dist_vec2(i,:)).^p)^2 ));
                        Denx(i) = (1 / ( (dist_vec2(i,j))^(p+2)  * sum(1./(dist_vec2(i,:)).^p)^2 ));               
                        Numy(i) = y_vec(i) * (1 / ( (dist_vec2(i,j))^(p+2)  * sum(1./(dist_vec2(i,:)).^p)^2 ));
                        Deny(i) = (1 / ( (dist_vec2(i,j))^(p+2)  * sum(1./(dist_vec2(i,:)).^p)^2 ));
                end
                XBS2(j) = sum(Numx)/sum(Denx);
                YBS2(j) = sum(Numy)/sum(Deny);
            end
        count = count + 1;
end
centroids(:,1) = XBS2;
centroids(:,2) = YBS2;