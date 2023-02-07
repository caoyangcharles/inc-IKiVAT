function [Pi_MST_1,length_partition_MST] = MST(cut,class,num_c)
[~,ind]=sort(cut,'descend');
IK_id_MST = 0;
IK_HL_MST = zeros(num_c-1,1);
ALL_label_MST = zeros(length(class),num_c-1);
for clusters = 2:num_c
    A = clusters-1;
    ind1=sort(ind(1:A));
    
    Pi_MST=zeros(length(class),1);
    Pi_MST(1:ind1(1)-1)=1;
    Pi_MST(ind1(end):end)=clusters;
    for k=2:clusters-1
        Pi_MST(ind1(k-1):ind1(k)-1)=k;
    end
    IK_id_MST = IK_id_MST + 1;
    ALL_label_MST(:,IK_id_MST) = Pi_MST;
    true_mem = class;
    mem = Pi_MST;
    AMI_MST = ami(true_mem,mem);
    IK_HL_MST(IK_id_MST,1) = AMI_MST;
end
[~,Par_MST] = max(IK_HL_MST);
Pi_MST_1 = ALL_label_MST(:,Par_MST);

length_partition_MST=zeros(1,length(unique(Pi_MST_1)));
for i=1:length(unique(Pi_MST_1))
    length_partition_MST(i)=length(find(Pi_MST_1==i));
end
end