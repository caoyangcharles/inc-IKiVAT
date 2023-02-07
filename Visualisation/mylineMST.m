function mylineMST(Pi_MST_1,length_partition_MST) 
     L1_MST = cumsum(length_partition_MST); 
    labnum1_MST = L1_MST(1);
    line([0.5,0.5],[0.5,labnum1_MST+0.5],'linestyle','--','color','b','LineWidth',2);
    line([0.5,labnum1_MST+0.5],[labnum1_MST+0.5,labnum1_MST+0.5],'linestyle','--','color','b','LineWidth',2);
    line([labnum1_MST+0.5,labnum1_MST+0.5],[0.5,labnum1_MST+0.5],'linestyle','--','color','b','LineWidth',2);
    line([0.5,labnum1_MST+0.5],[0.5,0.5],'linestyle','--','color','b','LineWidth',2);
    for i =1:length(unique(Pi_MST_1))-1
        labnum1_MST = L1_MST(i);
        labnum2_MST = L1_MST(i+1);
        line([labnum1_MST+0.5,labnum1_MST+0.5],[labnum1_MST+0.5,labnum2_MST+0.5],'linestyle','--','color','b','LineWidth',2);
        line([labnum1_MST+0.5,labnum2_MST+0.5],[labnum1_MST+0.5,labnum1_MST+0.5],'linestyle','--','color','b','LineWidth',2);
        line([labnum2_MST+0.5,labnum1_MST+0.5],[labnum2_MST+0.5,labnum2_MST+0.5],'linestyle','--','color','b','LineWidth',2);
        line([labnum2_MST+0.5,labnum2_MST+0.5],[labnum2_MST+0.5,labnum1_MST+0.5],'linestyle','--','color','b','LineWidth',2);
    end
  end