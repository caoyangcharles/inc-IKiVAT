function [X,AS,label_store] = Drawbar(I,class,colors)
% [X,AS,label_store] = Drawbar(smp_incsiVAT(I_incsiVAT),Labels);


%% 在inc-siVAT中 class指的是原本的数据标签，但是在这里传入的是整体的label，可是我们对此进行了采样，所以需要
%   使用采样后的样本的长度
[n,~] = size(I);
label_store = zeros(n,2);
B2 = 0; count1 =0;
for i = 1:n-1
    if class(I(i)) == class(I(i+1))
        B2 = B2+1;
    else
        B2 = B2+1;
        count1= count1+1;
        label_store(count1,1) = B2;
        label_store(count1,2) = class(I(i));
        hold on
        B2=0;
        continue;
    end
end
ls=sum(label_store);
if ls==n-1
    count1= count1+1;
    label_store(count1,1)=1;
    label_store(count1,2)=class(I(n));
else
    count1= count1+1;
    label_store(count1,1)=n-ls(1);
    label_store(count1,2)=class(I(n));
end
label_store(all(label_store == 0, 2),:) = [];
label_store = label_store';
[~,AS]=size(label_store);
X = zeros(AS,2)';
X(1,:) = label_store(1,:);
end