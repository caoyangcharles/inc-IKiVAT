function[RiV_incsiVAT,smp_incsiVAT,I_incsiVAT, cut_incsiVAT] =inc_siVAT(data_matrix_without_lables,cp, ns,times)

    total_no_of_points=size(data_matrix_without_lables,1); %%N
    dimension = size(data_matrix_without_lables,2); %%data dimension
    

    inc_sivat_window_size=500;  %%N_ch = InitN = 500
    no_of_iterations = floor(total_no_of_points/inc_sivat_window_size);
    
    data_till_now = data_matrix_without_lables(1:inc_sivat_window_size,:); %%X_curr
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%Steps 1 and 2, Apply siVAT (MMRS and iVAT) to the first or initial chunk
    [ RiV,rv,C,I,ri,cut,smp,m,Rp,nt_all,cp_nn_index,max_distance_Rp ] = siVAT(data_till_now, cp, ns, times); 
    
    cut=cut'; %%d - cut magnitudes
    cut(1)=[];

    rp_incsiVAT=Rp; %%R - distance matrix
    m_incsiVAT=m; %%M - Maximin points
    max_distance_Rp_incsiVAT=max_distance_Rp; %%Dmax

    RiV_incsiVAT=RiV;
    rv_incsiVAT=rv;
    C_incsiVAT=C; %%F = Connection Indices in MST
    I_incsiVAT=I; %% P - VAR Reordering
    ri_incsiVAT=ri; %%
    cut_incsiVAT=cut; %d - MST cut magnitudes
    smp_incsiVAT=smp; %%MMRS sample \tilde{S}
    %% 
    PerChunk_incsiVATTime=[];
    ElemstoRemove=[];
    ElemstoAdd=[];
    
    for i=2:no_of_iterations
    %%Step 3
    this_data=data_matrix_without_lables((i-1)*inc_sivat_window_size+1:i*inc_sivat_window_size,:);
    data_till_now_new=[data_till_now;this_data];
    
    [n,~]=size(data_till_now_new); %%N_curr
    nt_all_incsiVAT=zeros(1,cp); %%\mathcal{N}
    
    %%% Step 4
%     tic
    [m_incsiVAT,rp_incsiVAT,max_distance_Rp_incsiVAT] = incMMRS(data_till_now,this_data,cp,m_incsiVAT,rp_incsiVAT,max_distance_Rp_incsiVAT);
    [~,cp_nn_index_siVAT]=min(rp_incsiVAT,[],2); %% find  nearest neighbour to new cp points

%%% Step 5
            elemnts_to_remove=[];
            elements_to_add=[];

            tic
            for t=1:cp
                s = find(cp_nn_index_siVAT==t);
                nt = ceil(ns*length(s)/n);
                nt_all_incsiVAT(t)=nt;

                [tf, ~] = ismember(smp_incsiVAT, s);
                already_in_smp=smp_incsiVAT(tf);

                if(length(already_in_smp)>nt)
                    
                    rng(times)
                    ind=randperm(length(already_in_smp));
        
                    elemnts_to_remove=[elemnts_to_remove;already_in_smp(ind(1:length(already_in_smp)-nt))];
                else
                    if(length(already_in_smp)<nt)
                        [tf1, ~] = ismember(s,smp_incsiVAT);
                        not_in_smp=s(~tf1);
                        
                        rng(times);
                        ind=randperm(length(not_in_smp));
        
                        elements_to_add=[elements_to_add;not_in_smp(ind(1:nt-length(already_in_smp)))];
                    end
                end
            end
            
            
                        %%% Step 6  

            tic
            for j=1:length(elemnts_to_remove)
                point_to_remove=find(smp_incsiVAT==elemnts_to_remove(j));
                iVAT_point_to_remove_index=find(I_incsiVAT==point_to_remove);
                [rv_incsiVAT,C_incsiVAT,I_incsiVAT,ri_incsiVAT,cut_incsiVAT] = decVAT(rv_incsiVAT,C_incsiVAT,I_incsiVAT,ri_incsiVAT,cut_incsiVAT,point_to_remove);
                [RiV_incsiVAT] = deciVAT(rv_incsiVAT,RiV_incsiVAT,iVAT_point_to_remove_index);

                smp_incsiVAT_1=smp_incsiVAT;
                smp_incsiVAT_1(point_to_remove)=[];
                I_mapped=zeros(1,length(I_incsiVAT));
                for k=1:length(I_mapped)
                    idx=find(smp_incsiVAT_1==smp_incsiVAT(I_incsiVAT(k)));
                    I_mapped(k)=idx;
                end
                smp_incsiVAT=smp_incsiVAT_1;
                I_incsiVAT=I_mapped;
            end
            
             %%% Step 7

            for j=1:length(elements_to_add)
                distance_previous_points=distance2(data_till_now_new(elements_to_add(j),:),data_till_now_new(smp_incsiVAT(I_incsiVAT),:));
                [rv_incsiVAT,C_incsiVAT,I_incsiVAT,ri_incsiVAT,cut_incsiVAT,new_point_location] = incVAT(rv_incsiVAT,C_incsiVAT,I_incsiVAT,ri_incsiVAT,cut_incsiVAT,distance_previous_points);
                [RiV_incsiVAT] = inciVAT(rv_incsiVAT,RiV_incsiVAT,new_point_location);
                smp_incsiVAT=[smp_incsiVAT;elements_to_add(j)];
            end

            ElemstoRemove= [ElemstoRemove length(elemnts_to_remove)];
            ElemstoAdd= [ElemstoAdd length(elements_to_add)];


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Step 8
            data_till_now=data_till_now_new;
    end