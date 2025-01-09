function [Recommended_RN] = Meta_Distribution_function(system, reliability_threshod, MD_threshold)
%%
%System Parameters
Ps = system.Ps; %db2pow(23-30); %Transmit power (23dbm) in Watts for indoor nodes
Pm = system.Pm; %db2pow(40-30); %Transmit power (40dbm) in Watts for outdoor Macro cells
B = system.B; %100e6; %Total bandwidth per cell
No = system.No; %db2pow(-120 - 30)*B; %Noise density -120 dbm per Hz
fc = system.fc; %4e9; %Carrier frequency

%Propagation
K = system.K; %(3e8/(4*pi*fc))^2; %Path loss constant for indoor propagation (3GPP)
alpha = system.alpha; %3; %Path loss exponent for indoor propagation (3GPP) wthout wall material details

%Dimensions and Traffic Requirements
R_vec = system.R_vec; %linspace(5,15,30); %Radius of deployment area
H = system.H; %4; %Height of the floor
R_I = system.R_I; %20; %Interference range
N_u = system.N_u; %20; %Number of users
gamma = system.gamma; %db2pow(-10);

N_m = system.N_m; %3;
N_s_vector = system.N_s_vector; %1:20; %Number of indoor nodes
%Meta_dist = zeros(1, length(N_s_vector));


%Meta_dist = zeros(1, length(N_s_vector));

%Scenario Simulation
realizations = 10000;
    R = R_vec;
    lambda_u = N_u/(pi*R^2); %User density [users/m^2]
    lambda_a = 1e6; %User traffic [bits/s/user]
    sigma = lambda_a * lambda_u; %Traffic density [bits/s/m^2]
    Meta_dist = zeros(1, length(N_s_vector));
    flag = 0;

    for j = 1:length(N_s_vector)
        N_s = N_s_vector(j);
        for i = 1:realizations
            ang_user = 2*pi*rand;
            dist_user = R*sqrt(rand);
            x_user = dist_user*cos(ang_user);
            y_user = dist_user*sin(ang_user);

            dist_nodes_center = R*sqrt(rand(1, N_s));
            ang_nodes = 2*pi*rand(1, N_s);
            x_nodes = dist_nodes_center.*cos(ang_nodes);
            y_nodes = dist_nodes_center.*sin(ang_nodes);

            dist_macro_center = R + R_I*sqrt(rand(1, N_m));
            ang_macro = 2*pi*rand(1, N_m);
            x_macro = dist_macro_center.*cos(ang_macro);
            y_macro = dist_macro_center.*sin(ang_macro);

            dist_nodes_user = sqrt((x_nodes - x_user).^2 + (y_nodes- y_user).^2);
            dist_macro_user = sqrt((x_macro - x_user).^2 + (y_macro- y_user).^2);

            nearest_node_dist = min(dist_nodes_user);
            interfering_nodes_dist = dist_nodes_user(dist_nodes_user~=nearest_node_dist);

            Area_serving = pi*R^2/N_s;
            %gamma = 2^(lambda_a*lambda_u*Area_serving/B)-1;
            Term1 = exp(-(gamma*No)/(Ps * K * (H^2 + nearest_node_dist^2)^(-alpha/2)));
            Term2 = prod(1./(1 + gamma* (((H^2 + interfering_nodes_dist.^2).^(-alpha/2))/((H^2 + nearest_node_dist^2)^(-alpha/2)))));
            Term3 = prod(1./(1 + gamma* (((H^2 + dist_macro_user.^2).^(-alpha/2))/((H^2 + nearest_node_dist^2)^(-alpha/2)))));
            Conditional_success_interference = Term1*Term2*Term3; 
            Conditional_success_RSRP = Term1;
            if Conditional_success_RSRP > reliability_threshod
                Meta_dist(j) = Meta_dist(j) + 1/realizations;
            end
        end
        if Meta_dist(j) > MD_threshold
            Recommended_RN = N_s;
          break;
        end
    end