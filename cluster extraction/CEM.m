function [IP_CL,length_partition_IC] = CEM(RiV,class,num_c)
%% CEV  %%  Clustering based on images
result = zeros(length(RiV)-1,3);
cluster_set = zeros(length(RiV),length(RiV)-2);
% Calculate the distance between two adjacent rows
for i = 1:length(RiV)-1
    result(i,1) = pdist2(RiV(i,:),RiV(i+1,:));
    result(i,2)=i;
    result(i,3)=i+1;
end
merge_order = zeros(length(RiV)-1,3);
% There are a total of N categories at the beginning
clusters = zeros(length(RiV),1);
for j = 1:length(RiV)
    clusters(j) = j;
end
k =0;
% partition is the change in the number of clusters during the merging process
for partion = length(merge_order):-1:2
    [p,q] = min(result(:,1));
    merge_order(partion,1) = p;
    merge_order(partion,2) = result(q,2);
    merge_order(partion,3) = result(q,3);
    %% update the original cluster, the same category is represented by the same number
    if clusters(merge_order(partion,3)) ~= clusters(merge_order(partion,2))  
        clusters(clusters == clusters(merge_order(partion,3))) = clusters(merge_order(partion,2));
    end
    k = k+1;
    result(q,:) = [];
    %% Mark the number of clusters and store them in cluster_set
    QL = unique(clusters);
    QL(1) = [];
    cum = length(QL);
    Pi=zeros(length(RiV),1);
    
    Pi(1:QL(1))=1;
    Pi( QL(end):end)= cum+1;
    for L =2:cum
        Pi( QL(L-1): QL(L)-1)=L;
    end
    cluster_set(:,k) = Pi;
end
%%Search cluster
% label = zeros(length(RiV),1);
CL = zeros(num_c,1);
for C = 2: num_c
    lable  = cluster_set(:,length(RiV)-C);
    true_mem = class;
    mem = lable;
    fraction2 = ami(true_mem,mem);
    CL(C-1) = fraction2;
end
[~,W] = max(CL);  
IP_CL = cluster_set(:,length(RiV)-W-1);
%%   
length_partition_IC=zeros(1,length(unique(IP_CL))); 
for i=1:length(unique(IP_CL))
    length_partition_IC(i)=length(find(IP_CL==i)); 
end
end