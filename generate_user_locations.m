function  [x_vec,y_vec] = generate_user_locations(Xmax, Ymax, Xmin, Ymin, Nuser, optionx, optiony)


switch optionx
    case 1
        x_vec = Xmin + rand(1,Nuser)*(Xmax - Xmin); %Uniform
    case 2
        random2 =  exprnd(1, 1, Nuser);
        x_vec = Xmin + (random2/max(abs(random2)))*(Xmax - Xmin); %Exponentially Skewed
    case 3
        random3 = randn(1,Nuser);
        x_vec =(Xmax - Xmin)/2 *  (random3/max(abs(random3))); %Normal distributed
    case 4
        random41 =  exprnd(1, 1, ceil(Nuser/2));
        random42 = max(random41) - random41;
        x_vec = [Xmin + (random41/max(abs(random41)))*(Xmax - Xmin) Xmin + (random42/max(abs(random42)))*(Xmax - Xmin)]; %Biexponential
    otherwise
       disp('Invalid x input')
end

switch optiony
    case 1
        y_vec = Ymin + rand(1,Nuser)*(Ymax - Ymin); %Uniform
    case 2
        random2 =  exprnd(1, 1, Nuser);
        y_vec = Ymin + (random2/max(abs(random2)))*(Ymax - Ymin); %Exponentially Skewed
    case 3
        random3 = randn(1,Nuser);
        y_vec = (Ymax - Ymin)/2 *  (random3/max(abs(random3))); %Normal distributed
    case 4
        random41 =  exprnd(1, 1, ceil(Nuser/2));
        random42 = max(random41) - random41;
        y_vec = [Ymin + (random41/max(abs(random41)))*(Ymax - Ymin) Ymin + (random42/max(abs(random42)))*(Ymax - Ymin)];
    otherwise
       disp('Invalid x input')
end