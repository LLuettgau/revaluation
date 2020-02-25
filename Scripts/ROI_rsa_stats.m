    %script to run stats on rsa matrices for candyman
    
    clear all; close all; clc
    
    roi = 1; % 1 = L Hippocampus, 2 = R OFC
    
    cd '...'; %set path to /rsa_stats/
    
    %pre
        if roi == 1
            load('pre_pearson_sim_mat_L_Hippocampus_ROI.mat')
        elseif roi == 2
            load('pre_pearson_sim_mat_R_OFC_ROI.mat')
        end
    
    samples_pre = samples;
    
    ds_names_pre = {'S1A-US0'; 'S1A-US+'; ...
                    'S1B-US0'; 'S1B-US+'; ...
                    'S2A-US-'; 'S2A-US+'; ...
                    'S2B-US-'; 'S2B-US+'; ...
                    'S3A-US-'; 'S3A-US0'; ...
                    'S3B-US-'; 'S3B-US0'};

    figure;
    imagesc(samples); colorbar; hold all;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
   
   
    %post
    if roi == 1
       load('post_pearson_sim_mat_L_Hippocampus_ROI.mat')
    elseif roi == 2
       load('post_pearson_sim_mat_R_OFC_ROI.mat')
    end
    
    samples_post = samples;

    figure;
    imagesc(samples); colorbar; hold all;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
    

    %pre-post difference
    if roi == 1
       load('diff_pearson_sim_mat_L_Hippocampus_ROI.mat')
    elseif roi == 2
       load('diff_pearson_sim_mat_R_OFC_ROI.mat')
    end

    samples_prepost = samples;
    figure;
    imagesc(samples_prepost); colorbar; hold all;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
    
%     %pre-post change CS0A-US- and CS0B-US-, CS0A-US+ and CS0B-US+
%     (samples_prepost(3,1)+samples_prepost(4,2))/2
%     %2nd row in all_diff matrix  & 9th row in all_diff matrix 
    
    %t-test similarity changes against 0
        if roi == 1
            load('diff_pearson_all_corr_L_Hippocampus_ROI.mat')
        elseif roi == 2
            load('diff_pearson_all_corr_R_OFC_ROI.mat')
        end

%     %pre-post change CS+A-US- and CS+B-US-, CS+A-US0 and CS+B-US0
%     (samples_prepost(7,5)+samples_prepost(8,6))/2
%     %24th row in all_diff matrix  & 27th row in all_diff matrix 
    

    (samples_prepost(3,1)+samples_prepost(4,2))/2
    (samples_prepost(7,5)+samples_prepost(8,6))/2
    (samples_prepost(11,9)+samples_prepost(12,10))/2


    %pool similarity across all stimuli - 
    %CS-A/B --> 3,1 and 4,2 == 2nd and 13th row
    %CS0A/B --> 7,5 and 8,6 == 40th and 47th row
    %CS+A/B --> 11,9 and 12,10 == 62nd and 65th row
    
    diff_to_plot = [all_diff(2,:); all_diff(13,:); all_diff(40,:); all_diff(47,:); all_diff(62,:); all_diff(65,:)];
    pool_effect = mean(diff_to_plot(3:end,:));
    [h, p, ci, stats] = ttest(pool_effect,0, 'tail','left')
    [p, h, stats] = signrank(pool_effect,0)
    figure; hist(pool_effect)
    mean(pool_effect)
    std(pool_effect)
    
    pool_control = mean(diff_to_plot(1:2,:));
    [h, p, ci, stats] = ttest(pool_control,0, 'tail','right')
    [p, h, stats] = signrank(pool_control,0)
    figure; hist(pool_control)
    mean(pool_control)
    std(pool_control)
    
    %compare effect and control changes in rsa
    [h, p, ci, stats] = ttest(pool_effect,pool_control, 'tail','left')
    
    
    pool_csminus = mean(diff_to_plot(1:2,:));
    [h, p, ci, stats] = ttest(pool_csminus,0)
    [p, h, stats] = signrank(pool_csminus,0)
    mean(pool_csminus)
    std(pool_csminus)
    
    pool_csnull = mean(diff_to_plot(3:4,:));
    [h, p, ci, stats] = ttest(pool_csnull,0, 'tail','left')
    [p, h, stats] = signrank(pool_csnull,0)
    mean(pool_csnull)
    std(pool_csnull)
    
    pool_csplus = mean(diff_to_plot(5:6,:));
    [h, p, ci, stats] = ttest(pool_csplus,0, 'tail','left')
    [p, h, stats] = signrank(pool_csplus,0)
    mean(pool_csplus)
    std(pool_csplus) 
    
    %compare effect and control changes in rsa
    [h, p, ci, stats] = ttest(pool_csplus,pool_csminus, 'tail','left')
    [h, p, ci, stats] = ttest(pool_csnull,pool_csminus, 'tail','left')

    %CS-A-US0 and CS-B-US0
    [h, p, ci, stats] = ttest(diff_to_plot(1,:),0)
    mean(diff_to_plot(1,:))
    std(diff_to_plot(1,:)) 
    %CS-A-US+ and CS-B-US+
    [h, p, ci, stats] = ttest(diff_to_plot(2,:),0)
    mean(diff_to_plot(2,:))
    std(diff_to_plot(2,:)) 
    
    %CS0A-US- and CS0B-US-
    [h, p, ci, stats] = ttest(diff_to_plot(3,:),0, 'tail','left')
    mean(diff_to_plot(3,:))
    std(diff_to_plot(3,:)) 
    %CS0A-US+ and CS0B-US+
    [h, p, ci, stats] = ttest(diff_to_plot(4,:),0, 'tail','left')
    mean(diff_to_plot(4,:))
    std(diff_to_plot(4,:)) 


    %CS+A-US- and CS+B-US-
    [h, p, ci, stats] = ttest(diff_to_plot(5,:),0, 'tail','left')
    mean(diff_to_plot(5,:))
    std(diff_to_plot(5,:)) 
    %CS+A-US0 and CS+B-US0
    [h, p, ci, stats] = ttest(diff_to_plot(6,:),0, 'tail','left')
    mean(diff_to_plot(6,:))
    std(diff_to_plot(6,:)) 
    figure; hist(diff_to_plot(6,:))

    
    figure;
    subplot(131);  imagesc(samples_pre); colorbar; hold all;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
    title('pre')

    subplot(132);  imagesc(samples_post); colorbar;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
    title('post')
    
    subplot(133);  imagesc(samples_prepost); colorbar;
    set(gca,'XTick',1:12,'XTickLabel', ds_names_pre);
    set(gca,'YTick',1:12,'YTickLabel', ds_names_pre);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    set(gca,'TickDir','out')
    set(gca,'XTickLabelRotation', 45);
    title('pre-post change')
    
   