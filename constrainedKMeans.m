function [labels,centroids] = constrainedKMeans(X, K, tau, maxiter)
% Extract dataset size information
[n,D] = size(X);

% Parse inputs to algorithm
if nargin == 2
    tau = min(0.05*n, 0.5*n/K);
    maxiter = 100;
elseif nargin == 3
    maxiter=100;
end

% Assign initial clusters and centroids
labels = randi(K,n,1);
centroids = zeros(K,D);
for k = 1:K
    centroids(k,:) = mean(X(labels==k,:));
end

% Variable needed for while-loop.
iter = 1; % Used to ensure that we not exceed maxiter iterations. 

while iter<maxiter  
    
    % Below is our objective function. By default, this script uses 
    objectiveFunction = reshape(0.5*pdist2(X, centroids).^2, [n*K,1]); % Squared Euclidean distance between data points and centroids, vectorized

    % Ensure at most one cluster to be assigned to each point:  
    % A1 has exactly n*K nonzero entries in a (n,n*K) matrix.
    A1 = sparse(repmat((1:n)',1,K),  reshape(1:K*n, [n,K]), ones(n, K));
    b1 = ones(n,1); 

    % Ensure at least one cluster to be assigned to each point: 
    % A2 has exactly n*K nonzero entries in a (n,n*K) matrix.
    A2 = sparse(repmat((1:n)',1,K),  reshape(1:K*n, [n,K]), -ones(n, K));
    b2 = -ones(n,1); 

    % Enforce minimum cluster size constraint:
    % A2 has exactly n*K nonzero entries in a (K,n*K) matrix.
    A3 = sparse(reshape(repmat((1:K)', 1,n)', [n*K,1]),1:n*K, -ones(1,n*K)); 
    b3 = -tau*ones(K,1);

    % Run linear program to get new cluster assignments
    T = intlinprog(objectiveFunction, 1:n*K, [A1; A2; A3], [b1; b2; b3], [], [], zeros(n*K,1), ones(n*K,1));
    [~, labelsNew] = max(reshape(T,n,K), [], 2);

    % Assign new centroids: 
    centroidsNew = zeros(K,D);
    for k = 1:K
        centroidsNew(k,:) = mean(X(labelsNew == k,:));
    end

    % Compare centroids from prior iteration and current iteration. 
    % The following condition will be true if labels do not change across
    % an iteration. We stop in this case and output the current labels.
    if sum(diag(pdist2(centroids, centroidsNew)) == zeros(K,1)) == 2 % True if centroids do not change. 
        labels = labelsNew;
        centroids = centroidsNew;
        break
    else
        % Take current centroids and move to next iteration.
        centroids = centroidsNew;
        iter = iter+1;
    end
end
        


   