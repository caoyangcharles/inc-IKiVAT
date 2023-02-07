function mylineCER(IP_CL,length_partition_IC)
    L1_IC = cumsum(length_partition_IC); 

    labnum1_IC = L1_IC(1);
    line([0.5,0.5],[0.5,labnum1_IC+0.5],'linestyle','--','color','g','LineWidth',2);
    line([0.5,labnum1_IC+0.5],[labnum1_IC+0.5,labnum1_IC+0.5],'linestyle','--','color','g','LineWidth',2);
    line([labnum1_IC+0.5,labnum1_IC+0.5],[0.5,labnum1_IC+0.5],'linestyle','--','color','g','LineWidth',2);
    line([0.5,labnum1_IC+0.5],[0.5,0.5],'linestyle','--','color','g','LineWidth',2);
    for i =1:length(unique(IP_CL))-1
        labnum1_IC = L1_IC(i);
        labnum2_IC = L1_IC(i+1);
        line([labnum1_IC+0.5,labnum1_IC+0.5],[labnum1_IC+0.5,labnum2_IC+0.5],'linestyle','--','color','g','LineWidth',2);
        line([labnum1_IC+0.5,labnum2_IC+0.5],[labnum1_IC+0.5,labnum1_IC+0.5],'linestyle','--','color','g','LineWidth',2);
        line([labnum2_IC+0.5,labnum1_IC+0.5],[labnum2_IC+0.5,labnum2_IC+0.5],'linestyle','--','color','g','LineWidth',2);
        line([labnum2_IC+0.5,labnum2_IC+0.5],[labnum2_IC+0.5,labnum1_IC+0.5],'linestyle','--','color','g','LineWidth',2);
    end       
   end
