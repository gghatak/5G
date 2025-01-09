function [subgraphs] = k_flex_clustering(G, A, W, threshold)
subgraphs = G;
flag = 0;
Temp =[];
while flag == 0
    for i = 1: size(subgraphs, 1)
        if sum(W(subgraphs(i,subgraphs(i,:)~=0))) > threshold
            A_subgraph = A(subgraphs(i,subgraphs(i,:)~=0), subgraphs(i,subgraphs(i,:)~=0));
            [S, T] = MinCut([1], A_subgraph);
            First_indices = S(1,S(1,:)~=0);
            First_Group = subgraphs(i,First_indices);
            First_Group =[First_Group zeros(1, length(G) - length(First_Group))];
            Second_indices =  S(2,S(2,:)~=0);
            Second_Group = subgraphs(i,Second_indices);
            Second_Group =[Second_Group zeros(1, length(G) - length(Second_Group))];
            Temp = [Temp.',  First_Group.', Second_Group.'].';
        else
            Temp = [Temp.', subgraphs(i,:).'].';
        end
    end
    
    for i = 1:size(subgraphs,1)
        Wght(i) = sum(W(subgraphs(i, subgraphs(i,:)~=0)));
    end

    if sum(Wght > threshold) ==0
        flag = 1;
        break;
    end
    clear subgraphs
    subgraphs = Temp;
   Temp =[];
end
