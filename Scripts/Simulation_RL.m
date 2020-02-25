%Simulate choices based on best-fitting parameters

%==============================================================================
function []=Simulation_RL(version, mechanism)

% clear all; close all; clc

%version = 1; %indicate task version
pso = 0; %particle swarm optimization yes (1) or no (0)
%alphas = 2; %select number of learning rates, 1, 2, 3
inivals = 2; %set initial value for optimization: 1, 0, or 2 (based on grid search)
%mechanism = 1; %set value estimation algorithm, 1 = Rescorla-Wagner with 0, .5 and 1 as inputs, 2 = Rescorla-Wagner with 0/1 as inputs (depending on reward received or not),
%final value estimates are multiplied with individual rating
%of unconditioned stimulus during initial foodcue rating 
%3 = trial-by-trial scaling

num_simulations=10000;

% run version 1 (exp1 with parameters 1, 1)
% run version 2 (exp2 with parameters 2, 2)
% run version 3 (exp3 with parameters 3, 2)
% run version 4 (fmri with parameters 4, 2)

data_path = '/Users/lennart/behavioral_data/'; %set path to exp1, exp 2, exp3 or fmri here
res_dir =   '/Users/lennart/behavioral_data/sim_data/'; %specify folder to save structure with computational parameters

if version < 4
    params_dir= ['/Users/lennart/behavioral_data/params/candyman' num2str(version) filesep];
else
    params_dir= '/Users/lennart/behavioral_data/params/candyman_fmri/';
end
    
if version == 1
    rev_s = [3 5]; %exp1
    data_path = fullfile([data_path, 'exp1/']);
    datadir = data_path;
    vp = dir(fullfile([data_path, 'candymanLOG*']));
    load([params_dir 'cp_rw_3Alpha_1Tau_2_ini.mat']); %select best-fitting model
    load('/Users/lennart/behavioral_data/params/candyman1/identifier.mat'); %load identifier vector to only model included subjects
elseif version == 2
    rev_s = [1 3]; %exp2
    data_path = fullfile([data_path, 'exp2/']);
    datadir = data_path;
    vp = dir(fullfile([data_path, 'candyman2LOG*']));
    load([params_dir 'cp_asso_3Alpha_1Tau_2_ini.mat']); %select best-fitting model
    load('/Users/lennart/behavioral_data/params/candyman2/identifier.mat'); %load identifier vector to only model included subjects
elseif version == 3
    rev_s = [1 3 5]; %exp3
    data_path = fullfile([data_path, 'exp3/']);
    datadir = data_path;
    vp = dir(fullfile([data_path, 'candyman3LOG*']));
    load([params_dir 'cp_asso_3Alpha_1Tau_2_ini.mat']); %select best-fitting model
    load('/Users/lennart/behavioral_data/params/candyman3/identifier.mat'); %load identifier vector to only model included subjects
elseif version == 4
    rev_s = [3 5]; %fmri
    data_path = fullfile([data_path, 'fmri/']);
    datadir = data_path;
    vp = dir(fullfile([data_path, 'candyman_fmriLOG_9*']));
    load([params_dir 'cp_asso_3Alpha_1Tau_2_ini.mat']); %select best-fitting model
    load('/Users/lennart/behavioral_data/params/candyman_fmri/identifier.mat'); %load identifier vector to only model included subjects
end

vp = vp(identifier);
vp = {vp.name};

if version == 4
   predatadir = data_path;
   predata = dir(fullfile([data_path,'candyman_fmri_preLOG_9*']));
   predata = {predata.name};
   predata = predata(identifier);
end

sim_choices = NaN(120,num_simulations,numel(vp));
decision_trials = NaN(120,2,numel(vp));

for iteration=1:numel(vp)
    disp('Current iteration is:'); disp(iteration); drawnow
    
    
    %% preparations
    clear trials
    clear r
    clear r_new
    clear r_rev
    clear data
    
    tic %timing for each participant
    data = importdata([datadir, vp{iteration}]);
    
    name = data.infoP;
    
    fit_params=cp.params(iteration,:);
    
    if version == 4
        FOC = importdata([predatadir, predata{iteration}]);
        FOC = FOC.FOC;
        data.FOC = FOC;
    end
    
    %delete NaNs from devaluation
    if version == 3
        posnan = find(data.devaluation(:,2) == 1 | data.devaluation(:,2) == 2 | data.devaluation(:,2) == 3);
    else
        posnan = find(data.devaluation(:,2) == 1 | data.devaluation(:,2) == 2);
    end
    
    devaluation = data.devaluation;
    devaluation = devaluation(posnan,:);
    devalpos = find(isnan(devaluation(:,12)));
    devaluation(devalpos,:) = [];
    
    %define choice R/L, where the CS associated with the higher valued US was
    %chosen
    highvalchoiceL = length(find(devaluation(:,2) > devaluation(:,3) & devaluation(:,12) == 1));
    highvalchoiceR = length(find(devaluation(:,3) > devaluation(:,2) & devaluation(:,12) == 2));
    
    %percent choices of higher valued CS - devaluation --> manipulation check!
    percentchoicehighval = (highvalchoiceL + highvalchoiceR)/length(devaluation);
    results(iteration,1) = percentchoicehighval;
    
    if results(iteration,1) >= .5 && sum(data.FOC(:,4) > 0) > length(data.FOC(:,4))-(length(data.FOC(:,4))*.90)
        
        %extract revaluation choice data
        revaluation = data.devaluation;
        if version == 1 || version == 2 || version == 4
            revaluation_pos = revaluation(revaluation(:,2) == 1 | revaluation(:,2) == 2);
            revaluation = revaluation(revaluation_pos,:);
        elseif version == 3
            revaluation_pos = revaluation(revaluation(:,2) == 1 | revaluation(:,2) == 2 | revaluation(:,2) == 3);
            revaluation = revaluation(revaluation_pos,:);
        end
        revaluation = revaluation(~isnan(revaluation(:,12)),:);
        
        %extract forced choice data to generate LLE
        forcedchoices = data.forcedchoicekanjis(:,12);
        choice_trials = data.forcedchoicekanjis(:,2:3);
        choice_trials = choice_trials(~isnan(forcedchoices),:);
        forcedchoices = forcedchoices(~isnan(forcedchoices));
        
        %find positions of different S (1-6) during learning phase
        for k = 1:6
            kanjis_pos = find(data.forcedchoicekanjis(:,2)==k);
            kanjis(k) = data.forcedchoicekanjis(kanjis_pos(1),4);
            all_kanjis = data.FOC(:,5);
            trials(:,k) = find(all_kanjis == kanjis(k));
        end
        %% learning phase - gives an values estimate for all 6 CS (alpha1)
        data = data.FOC;
        
        r = data(:,8); %received inputs (rewards)
        
        %recode rewards - 1 = association with most appetitive reward, 0.5 = intermediate reward, 0 = least appetitive reward
        unir = unique(r);
        r(r==(unir(1))) = 0;
        r(r==(unir(2))) = 0.5;
        r(r==(unir(3))) = 1;
        r(data(:,14)==0) = NaN; %if no reward was presented (probabilistic feedback), value should not be updated
        
        sub_r = data(:,8); %received inputs (subjective rewards)
        
        %assess initial ratings
        uni_sub_r = unique(sub_r)';
        
        %standardize rating (0-1)
        u_sub_r(1:2) = uni_sub_r(1)/100;
        u_sub_r(3:4) = uni_sub_r(2)/100;
        u_sub_r(5:6) = uni_sub_r(3)/100;
        
        if mechanism == 2 || mechanism == 3
            r = data(:,14); %reward presented or not?
        end
        
        for k = 1:6
            r_new(:,k) = r(trials(:,k));
        end
        r = r_new;
        
        
        %preparations for revaluation phase
        choices = revaluation(:,12); %choices during revaluation, 1 left, 2 right
        
        if version == 1 || version == 2 || version == 4
            opt_order = revaluation(:,2); %during revaluation, 1 left/right, 2 right/left
            for z = 1:size(choices,1)
                if choices(z) == 1 && opt_order(z) == 1 || choices(z) == 2 && opt_order(z) == 2
                    r_rev(z,:) = [1 -1]; %update value of chosen option by +1, update unchosen option by -1
                elseif choices(z) == 2 && opt_order(z) == 1 || choices(z) == 1 && opt_order(z) == 2
                    r_rev(z,:) = [-1 1]; %update value of chosen option by +1, update unchosen option by -1
                end
            end
        else
            opt_order = revaluation(:,2:3); %during revaluation, 1 left/right, 2 right/left
            for z = 1:size(choices,1)
                if sum([choices(z) == 1 opt_order(z,:) == [1 2]]) == 3 || sum([choices(z) == 2  opt_order(z,:) == [2 1]]) == 3
                    r_rev(z,:) = [1 -1 NaN]; %update value of chosen option by +1, update unchosen option by -1
                elseif sum([choices(z) == 2 opt_order(z,:) == [1 2]]) == 3 || sum([choices(z) == 1 opt_order(z,:) == [2 1]]) == 3
                    r_rev(z,:) = [-1 1 NaN]; %update value of chosen option by +1, update unchosen option by -1
                elseif sum([choices(z) == 1 opt_order(z,:) == [2 3]]) == 3 || sum([choices(z) == 2 opt_order(z,:) == [3 2]]) == 3
                    r_rev(z,:) = [NaN 1 -1]; %update value of chosen option by +1, update unchosen option by -1
                elseif sum([choices(z) == 2 opt_order(z,:) == [2 3]]) == 3 || sum([choices(z) == 1 opt_order(z,:) == [3 2]]) == 3
                    r_rev(z,:) = [NaN -1 1]; %update value of chosen option by +1, update unchosen option by -1
                end
            end
        end
        
        for ridx=1:num_simulations
            [sim_choices(:,ridx,iteration)]=make_choices_candyman(fit_params,version,mechanism, u_sub_r, r, rev_s, r_rev, choices, choice_trials, forcedchoices);
        end
        
        %% compare simulated choices with forced choices        
        %loop over all CS for all binary combinations of choices and compute choice probability
        
        %order is 1 = CS-A, 2 = CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B

        for k = 1:num_simulations
            for m = 1:6 %order is 1 = CS-A, 2 = CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
                [posL] = find(choice_trials(:,1) == m);
                [posR] = find(choice_trials(:,2) == m);
                CS(k,m) = (sum(sim_choices(posL,k,iteration) == 1) + sum(sim_choices(posR,k,iteration) == 2)) / (length(posL) + length(posR));
            end
        end
        
        %order is CS-A, CS-B, CS0A, CS0B, CS+A, CS+B
        sim_results(iteration,1) = name;
        sim_results(iteration,2:7) = mean(CS);
        
    end    
end

    %save simulation results
    cd(res_dir)
	if version == 1
        save([pwd,filesep,'exp1_simulated_3Alpha_1Tau_RW.mat'],'sim_choices','CS','sim_results');
    elseif version == 2
        save([pwd,filesep,'exp2_simulated_3Alpha_1Tau_Asso.mat'],'sim_choices','CS','sim_results');
    elseif version == 3
        save([pwd,filesep,'exp3_simulated_3Alpha_1Tau_Asso.mat'],'sim_choices','CS','sim_results');
    elseif version == 4
        save([pwd,filesep,'fmri_simulated_3Alpha_1Tau_Asso.mat'],'sim_choices','CS','sim_results');
    end


res_to_plot = sim_results(:,2:end);

median(res_to_plot)
max(res_to_plot)
min(res_to_plot)

% repeated measures ANOVA on overall choice probabilities
data = res_to_plot;
varNames = {'S1A','S1B','S2A','S2B','S3A','S3B'};

t = array2table(data,'VariableNames',varNames);

factorNames = {'valence','stimulus'};
within = table({'low'; 'low';'med';'med';'hi';'hi'},{'A';'B';'A';'B';'A';'B'},'VariableNames',factorNames);


rm = fitrm(t,'S1A-S3B~1','WithinDesign',within); 


[ranovatbl] = ranova(rm, 'WithinModel', 'valence*stimulus');

disp('overall choice probability rmANOVA')
ranovatbl

% Wilcoxon signed-rank test on overall choice probabilities
%CS1A = CS-A vs. CS1B = CS-B

disp('overall choice probability CS-A vs. CS-B')
[P,H,STATS] = signrank(res_to_plot(:,1), res_to_plot(:,2))

%CS2A = CS0A vs. CS2B = CS0B
disp('overall choice probability CS0A vs. CS0B')
[P,H,STATS] = signrank(res_to_plot(:,3), res_to_plot(:,4))


%CS3A = CS+A vs. CS3B = CS+B
disp('overall choice probability CS+A vs. CS+B')
[P,H,STATS] = signrank(res_to_plot(:,5), res_to_plot(:,6))




%% plot simulated behavioral results
%used for Figure S3

cols.k = [0 0 0];
cols.b = [0   15 175];
cols.y = [255 211 0];

color_mat = [cols.k; cols.k; cols.y; cols.k; cols.b; cols.k];

cols.k = [0 0 0];
cols.b = [0 .058 .686];
cols.y = [1  .828 0];
cols.grey = [0.7843 0.7843 0.7843];
cols.dgrey = [0.1922 0.2000 0.2078];

positions = [1 2 4 5 7 8];   
pos = [positions-.15; ...
       positions+.15];

if version == 1
    %exp1
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0); hold all; 
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k); hold all;
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    LabelsCS ={'CS^{-}_{A}', 'CS^{-}_{B}', 'CS^{0}_{A}', 'CS^{0}_{B}', 'CS^{+}_{A}', 'CS^{+}_{B}'};
    ylim([0 1.05]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.25:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',positions)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    %set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    %ylabel('Choice Probability', 'FontSize',60,'FontType', 'Arial','Color','k') 
      % prepend a color for each tick label
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new);
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    %set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew);
keyboard
    
elseif version == 2
    %exp2
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0); hold all;
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y); hold all;
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    LabelsCS ={'CS^{-}_{A}', 'CS^{-}_{B}', 'CS^{0}_{A}', 'CS^{0}_{B}', 'CS^{+}_{A}', 'CS^{+}_{B}'};
    ylim([0 1.05]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.25:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',positions)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    %set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    %ylabel('Choice Probability', 'FontSize',60,'FontType', 'Arial','Color','k') 
      % prepend a color for each tick label
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new);
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew);
keyboard
elseif version == 3
    %exp3
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0);hold all;
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y); hold all;
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.b);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    LabelsCS ={'CS^{-}_{A}', 'CS^{-}_{B}', 'CS^{0}_{A}', 'CS^{0}_{B}', 'CS^{+}_{A}', 'CS^{+}_{B}'};
    ylim([0 1.05]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.25:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',positions)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    %set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    %ylabel('Choice Probability', 'FontSize',60,'FontType', 'Arial','Color','k') 
      % prepend a color for each tick label
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new);
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew);

    keyboard


elseif version == 4
    %fmri
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0);hold all;
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k); hold all;
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k);
    LabelsCS ={'CS^{-}_{A}', 'CS^{-}_{B}', 'CS^{0}_{A}', 'CS^{0}_{B}', 'CS^{+}_{A}', 'CS^{+}_{B}'};
    ylim([0 1.05]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.25:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',positions)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    %set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    %ylabel('Choice Probability', 'FontSize',60,'FontType', 'Arial','Color','k') 
      % prepend a color for each tick label
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new);
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew);
keyboard

end
