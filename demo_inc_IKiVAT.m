clear
close all

%%%%%%:Load Input (simulated) data stream %%%%%%%%%%%%%%%%%
%% %%%%%  SCData1.mat  %%%%%
load('SCData1.mat');
load('colors.mat');

times = 6;  
cp = 20;
ns = 150;
ba = 10;

%%  Reorder and normalize data
rng(times);
p = randperm(size(data_matrix_without_lables,1));    
data_matrix_without_lables = data_matrix_without_lables(p,:);
Labels = Labels(p,:);

data_matrix_without_lables = (data_matrix_without_lables - min(data_matrix_without_lables)).*((max(data_matrix_without_lables) - min(data_matrix_without_lables)).^-1);
data_matrix_without_lables(isnan(data_matrix_without_lables)) = 0.5;

%% Map to Isolation kernel space
tic
psi = 32;
rng(times);
[data_matrix_without_lables] = IKspace(data_matrix_without_lables(1:min(2000,size(Labels)),:),data_matrix_without_lables, psi, 200);
time_IKspace=toc;
%%  %%%%% inc-siVAT %%%%%

[RiV_incsiVAT,smp_incsiVAT,I_incsiVAT, cut_incsiVAT] = inc_siVAT(data_matrix_without_lables,cp, ns,times);
num_c = ceil(sqrt(length(RiV_incsiVAT)));


%%%%%%%%%%%%%%%%%%%%%%%%%%% Use CER methods to extract clusters from the RDI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[IP_IC,length_partition_IC] = CEM(RiV_incsiVAT,Labels(smp_incsiVAT(I_incsiVAT)),num_c);
[Pi_MST_1,length_partition_MST] = MST(cut_incsiVAT,Labels(smp_incsiVAT(I_incsiVAT)),num_c);

true_mem = Labels(smp_incsiVAT(I_incsiVAT));
score_CEM = ami(true_mem, IP_IC)     % see Isolation kernel CEM in Table 1 of paper
score_MST = ami(true_mem, Pi_MST_1)  % see Isolation kernel MST in Table 1 of paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw RDI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
figure
b = axes('Units', 'Normalized', 'Position', [ 0.1 0.1 0.7 0.7]);
imagesc(RiV_incsiVAT); colormap(gray);  axis image;axis off;


mylineCER(IP_IC,length_partition_IC); % draw cluster boundaries
mylineMST(Pi_MST_1,length_partition_MST);% draw cluster boundaries
%%%%%%%% Add a bar to show ground truth cluster labels. %%%%%%%
[X,AS,label_store] = Drawbar(smp_incsiVAT(I_incsiVAT),Labels);%%% 
axes('Units', 'Normalized', 'Position', [0.7 0.1 0.1 0.7]);  %%% 
b = bar(X,ba,'stack','EdgeColor','none');set(gca,'ydir','reverse');
axis image;axis off;

for n=1:AS
    b(n).FaceColor = colors(label_store(2,n),:);
end 