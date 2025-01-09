function  [labels,centroids] = WeightedKHarmonicMeans(X, K, maxiter)

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
q = 1;
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

                        T1 = sum((q./dist_vec2(i,:).^((p+q)))) * dist_vec2(i,j).^(-(q+2));
                        T2 = ((p+q)./dist_vec2(i,j).^(p+q+2)) * sum(dist_vec2(i,:).^(-(p+q)));
                        T3 = ( sum(dist_vec2(i,:).^(-q)) * sum(dist_vec2(i,:).^(-(p+q))) )^2;

                        Numx(i) = x_vec(i) * (T1 + T2)/T3 ;                 
                        Denx(i) =  (T1 + T2)/T3;               
                        Numy(i) = y_vec(i) *  (T1 + T2)/T3;
                        Deny(i) =  (T1 + T2)/T3;
                end
                XBS2(j) = sum(Numx)/sum(Denx);
                YBS2(j) = sum(Numy)/sum(Deny);
                if (isnan(XBS2(j)) || isnan(YBS2(j)))
                    flag = 0;
                end
            end
        count = count + 1;
end
centroids(:,1) = XBS2;
centroids(:,2) = YBS2;