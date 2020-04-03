%Analysis script for Experiment 5
%(Luettgau, Tempelmann, Kaiser, Jocham)

clc; clear all; close all;

plots  = 1;

exclude_outliers = 1; if exclude_outliers; choice_crit = .50; end %0 = all, 1 = no outliers, 2 = only outliers 

addpath(genpath('...')) %set path to Effect size toolbox (needs to be downloaded at: https://github.com/hhentschke/measures-of-effect-size-toolbox)

cd '...' %set path to exp5

data = dir('candyman_controlLOG_4*.mat');

format compact 


for i = 1:length(data)
load(data(i).name)

name = candyman_control.infoP;


results(i,1) = name;

    
    %%range of selected kanjis
    range_kanjis(i,1) = min(candyman_control.ratingkanjis(2,1:4));
    range_kanjis(i,2) = max(candyman_control.ratingkanjis(2,1:4));
    
    %how many button presses during learning phase (cover task) 
    if sum(candyman_control.FOC(:,3) == 0) > 160 -(160*.10)
       results(i,146) = 1;
    end
       results(i,145) = sum(candyman_control.FOC(:,3) == 0);


%% Manipulation check


%delete NaNs from devaluation
posnan = find(candyman_control.devaluation(:,2) == 1 | candyman_control.devaluation(:,3) == 1 ...
              | candyman_control.devaluation(:,2) == 2 | candyman_control.devaluation(:,3) == 2 ...
              | candyman_control.devaluation(:,2) == 3 | candyman_control.devaluation(:,3) == 3 ...
              | candyman_control.devaluation(:,2) == 4 | candyman_control.devaluation(:,3) == 4);
devaluation = candyman_control.devaluation;
devaluation = devaluation(posnan,:);
devalpos = find(isnan(devaluation(:,12)));
devaluation(devalpos,:) = [];

%define choice R/L, where the CS associated with the higher valued US was
%chosen
highvalchoiceL = length(find(devaluation(:,2) > devaluation(:,3) & devaluation(:,12) == 1));
highvalchoiceR = length(find(devaluation(:,2) < devaluation(:,3) & devaluation(:,12) == 2));

%percent choices of higher valued CS - devaluation --> manipulation check!
percentchoicehighval = (highvalchoiceL + highvalchoiceR)/length(devaluation);
results(i,2) = percentchoicehighval;


%percent choices for CS0A vs CS+A and CS0B vs CS+B

%calculate choice prob. between CS0A vs CS+A (80%)
CS3AvsCS1AchoiceL = length(find(devaluation(:,2) == 3 & devaluation(:,3) == 1 & devaluation(:,12) == 1));
CS3AvsCS1AchoiceR = length(find(devaluation(:,2) == 1 & devaluation(:,3) == 3 & devaluation(:,12) == 2));
percentchoiceCS3AvsCS1A = (CS3AvsCS1AchoiceL+CS3AvsCS1AchoiceR)/(length(find(devaluation(:,2) == 1 & devaluation(:,3) == 3)) + length(find(devaluation(:,2) == 3 & devaluation(:,3) == 1)));
results(i,143) = percentchoiceCS3AvsCS1A;
        
%calculate choice prob. between CS0B vs CS+B (20 %)
CS3BvsCS1BchoiceL = length(find(devaluation(:,2) == 4 & devaluation(:,3) == 2 & devaluation(:,12) == 1));
CS3BvsCS1BchoiceR = length(find(devaluation(:,2) == 2 & devaluation(:,3) == 4 & devaluation(:,12) == 2));
percentchoiceCS3BvsCS1B = (CS3BvsCS1BchoiceL + CS3BvsCS1BchoiceR)/(length(find(devaluation(:,2) == 4 & devaluation(:,3) == 2)) + length(find(devaluation(:,2) == 2 & devaluation(:,3) == 4)));
results(i,144) = percentchoiceCS3BvsCS1B;


%RT for higher associated CS - devaluation 
highvalL = find(devaluation(:,2) == 1 & devaluation(:,3) == 3);
highvalR = find(devaluation(:,2) == 3 & devaluation(:,3) == 1);
highvalRT = [devaluation(highvalL,13);
             devaluation(highvalR,13)];
         
%mean RT for higher associated CS - devaluation 
medianRThighval = nanmedian(devaluation([highvalL; highvalR], 13)); 
results(i,3) =  medianRThighval;

%RT for lower associated CS  - devaluation 
lowvalL = find(devaluation(:,2) == 2 & devaluation(:,3) == 4);
lowvalR = find(devaluation(:,2) == 4 & devaluation(:,3) == 2);
lowvalRT = [devaluation(lowvalL,13);
            devaluation(lowvalR,13)];
        
%mean RT for lower associated CS - devaluation 
medianRTlowval = nanmedian(devaluation([lowvalL; lowvalR], 13)); 
results(i,4) =  medianRTlowval;

results(:,147) = results(:,1);


%% Forced Choice Probe

%delete NaNs from forced choice responses (kanjis)
real_decisions = find(candyman_control.forcedchoicekanjis(:,2) == 1 | candyman_control.forcedchoicekanjis(:,3) == 1 |...
                      candyman_control.forcedchoicekanjis(:,2) == 2 | candyman_control.forcedchoicekanjis(:,3) == 2 |...
                      candyman_control.forcedchoicekanjis(:,2) == 3 | candyman_control.forcedchoicekanjis(:,3) == 3 |...
                      candyman_control.forcedchoicekanjis(:,2) == 4 | candyman_control.forcedchoicekanjis(:,3) == 4);
                  
forcedchoicekanjis = candyman_control.forcedchoicekanjis;
forcedchoicekanjis = forcedchoicekanjis(real_decisions,:);
posnankanjis = find(isnan(forcedchoicekanjis(:,12)));

forcedchoicekanjis(posnankanjis,:) = [];


%% loop over all stimuli
for m = 1:4
    [posL] = find(forcedchoicekanjis(:,2) == m);
    [posR] = find(forcedchoicekanjis(:,3) == m);

    CS(i,m) = (sum(forcedchoicekanjis(posL,12) == 1) + sum(forcedchoicekanjis(posR,12) == 2)) / (length(posL) + length(posR));
end

%% dumb coding
%percent CS1A chosen
[pos_cS1AL] = find(forcedchoicekanjis(:,2) == 1);
[pos_cS1AR] = find(forcedchoicekanjis(:,3) == 1);
[cS1AL] = forcedchoicekanjis(pos_cS1AL,4);
results(i,6) = cS1AL(1);

CS1A = sum(forcedchoicekanjis(pos_cS1AL,12) == 1) + sum(forcedchoicekanjis(pos_cS1AR,12) == 2);
percentchoiceCS1A = CS1A/(length(pos_cS1AL) + length(pos_cS1AR));%length(forcedchoicekanjis(pos_cS1AAL,12) == 1 & forcedchoicekanjis(pos_cS1AAR,12) == 2);
results(i,10) = percentchoiceCS1A;

%percent CS1B chosen
[pos_cS1BL] = find(forcedchoicekanjis(:,2) == 2);
[pos_cS1BR] = find(forcedchoicekanjis(:,3) == 2);
[cS1BL] = forcedchoicekanjis(pos_cS1BL,4);
results(i,7) = cS1BL(1);

CS1B = sum(forcedchoicekanjis(pos_cS1BL,12) == 1) + sum(forcedchoicekanjis(pos_cS1BR,12) == 2);
percentchoiceCS1B = CS1B/(length(pos_cS1BL) + length(pos_cS1BR)); %length(forcedchoicekanjis(pos_cS1BL,12) == 1 & forcedchoicekanjis(pos_cS1BR,12) == 2);
results(i,11) = percentchoiceCS1B;


%percent CS3A chosen
[pos_cS3AL] = find(forcedchoicekanjis(:,2) == 3);
[pos_cS3AR] = find(forcedchoicekanjis(:,3) == 3);
[cS3AL] = forcedchoicekanjis(pos_cS3AL,4);
results(i,8) = cS3AL(1);

CS3A = sum(forcedchoicekanjis(pos_cS3AL,12) == 1) + sum(forcedchoicekanjis(pos_cS3AR,12) == 2);
percentchoiceCS3A = CS3A/(length(pos_cS3AL) + length(pos_cS3AR)); % forcedchoicekanjis(pos_cS3NDL,12) == 1 & forcedchoicekanjis(pos_cS3NDR,12) == 2);
results(i,12) = percentchoiceCS3A;

%percent CS3B chosen
[pos_cS3BL] = find(forcedchoicekanjis(:,2) == 4);
[pos_cS3BR] = find(forcedchoicekanjis(:,3) == 4);
[cS3BL] = forcedchoicekanjis(pos_cS3BL,4);
results(i,9) = cS3BL(1);

CS3B = sum(forcedchoicekanjis(pos_cS3BL,12) == 1) + sum(forcedchoicekanjis(pos_cS3BR,12) == 2);
percentchoiceCS3B = CS3B/(length(pos_cS3BL) + length(pos_cS3BR)); %length(forcedchoicekanjis(pos_cS3BL,12) == 1 & forcedchoicekanjis(pos_cS3BR,12) == 2);
results(i,13) = percentchoiceCS3B;

%% Reaction times during decision probe

%RT for CS+ decisions
highvalL = find(forcedchoicekanjis(:,2) == 3 & forcedchoicekanjis(:,3) == 4);
highvalR = find(forcedchoicekanjis(:,2) == 4 & forcedchoicekanjis(:,3) == 3);
highvalRT = [forcedchoicekanjis(highvalL,13);
             forcedchoicekanjis(highvalR,13)];
         
%mean RT for higher valued CS chosen - devaluation 
medianRThighval = nanmedian(forcedchoicekanjis([highvalL; highvalR], 13)); 
results(i,15) =  medianRThighval;

%RT for lower valued CS chosen - devaluation 
lowvalL = find(forcedchoicekanjis(:,2) == 1 & forcedchoicekanjis(:,3) == 2);
lowvalR = find(forcedchoicekanjis(:,2) == 2 & forcedchoicekanjis(:,3) == 1);
lowvalRT = [devaluation(lowvalL,13);
            devaluation(lowvalR,13)];
        
%mean RT for lower valued CS chosen - devaluation 
medianRTlowval = nanmedian(devaluation([lowvalL; lowvalR], 13)); 
results(i,16) =  medianRTlowval;


%% Difference score pre- to post-exp rating of CS

for j = 1:4
    prepostkanjis = find(candyman_control.ratingkanjis(5,j) == candyman_control.postratingkanjis(5,:));
    value = find(candyman_control.ratingkanjis(5,j) == candyman_control.FOC(:,5));
    valuekanji = candyman_control.FOC(value(1),8);
    order_kanjis = candyman_control.FOC(value(1),26);
    
    diffkanji(1,j) = candyman_control.ratingkanjis(2,j);
    diffkanji(2,j) = candyman_control.postratingkanjis(2,prepostkanjis);
    diffkanji(3,j) = candyman_control.postratingkanjis(2,prepostkanjis) - candyman_control.ratingkanjis(2,j);
    diffkanji(4,j) = valuekanji;
    diffkanji(5,j) = candyman_control.FOC(value(1),5);
    diffkanji(6,j) = order_kanjis;
    
end
diffkanji = diffkanji';
diffkanji = sortrows(diffkanji,6);
diffkanji = diffkanji';
results(i,48:51) = diffkanji(1,:);

results(i,52:55) = diffkanji(2,:);

results(i,56:59) = diffkanji(3,:);

results(i,60:63) = diffkanji(4,:);

j = [];

%% Difference score pre- to post-exp rating of US

if isfield(candyman_control,'postratingfoodcues')

for j = 1:2
    
    val =  unique(candyman_control.FOC(:,7));
    prefood = find(candyman_control.ratingfoodcues(5,:) == val(j));
    postfood = find(candyman_control.postratingfoodcues(5,:) == val(j));
   
    difffood(1,j) = val(j);
    difffood(2,j) = candyman_control.ratingfoodcues(2,prefood);
    difffood(3,j) = candyman_control.postratingfoodcues(2,postfood);
    difffood(4,j) = candyman_control.postratingfoodcues(2,postfood) - candyman_control.ratingfoodcues(2,prefood);

end
difffood = difffood';
difffood = sortrows(difffood,2);
difffood = difffood';
results(i,64:65) = difffood(1,:);

results(i,66:67) = difffood(2,:);

results(i,68:69) = difffood(3,:);

results(i,70:71) = difffood(4,:);

else
results(i,64:71) = NaN;

end

if exclude_outliers == 1
    if results(i,2) < choice_crit
       results(i,2:end) = NaN;
       range_kanjis(i,:) = NaN;
    end
elseif exclude_outliers == 2
    if results(i,2) > .5
       results(i,2:end) = NaN;
    end
end


end

if exclude_outliers == 1
    results = results(~isnan(results(:,2)),:);
elseif exclude_outliers == 2
    results = results(~isnan(results(:,2)),:);   
end


%% CS/US subjective value (pre-rating)   

%Pre-ratings of US
pre_food = results(:,66:67);
nanmean(pre_food)
nanstd(pre_food)

%Pre-ratings of CS
disp('range of pre Kanjis ratings')
nanmean(range_kanjis)
nanmean(range_kanjis(:,2) - range_kanjis(:,1)) 

% repeated measures ANOVA on CS subjective value (pre-rating) 
pre_kanji = results(:,48:51);
data = pre_kanji;
varNames = {'S1','S2','S3','S4'};

t = array2table(data,'VariableNames',varNames);

factorNames = {'stimulus'};
within = table({'one'; 'two';'three';'four'},'VariableNames',factorNames);


rm = fitrm(t,'S1-S4~1','WithinDesign',within); 


[ranovatbl] = ranova(rm, 'WithinModel', 'stimulus');

disp('rmANOVA on kanjis ratings')
ranovatbl


%calculate partial eta square 
disp('overall rmANOVA effect size of ME stimulus pre')
data_rm = [data(:,1) data(:,2) data(:,3) data(:,4)];

mes1way(data_rm,'partialeta2','isDep',1)


%% overall choice probabilities for CS during choice probe phase
%terminology is CS1A = CS080, CS1B = CS020, CS3A = CS+80, CS3B = CS+20

cS1A = results(:,10);
cS1B = results(:,11);
cS3A = results(:,12);
cS3B = results(:,13);

res_to_plot = [cS1A cS1B cS3A cS3B];

mean(res_to_plot)
std(res_to_plot)

median(res_to_plot)
max(res_to_plot)
min(res_to_plot)

%compare revaluation choice probabilities
res_deval = [results(:,2) results(:,143:144)];
[p,h,STATS] = signrank(res_deval(:,2), res_deval(:,3))
mes(res_deval(:,2), res_deval(:,3),'U3')

median(results(:,143:144))
max(results(:,143:144))
min(results(:,143:144))

disp('all S')
[h,p,ci,stats] = ttest(res_to_plot, 0.5)

disp('S1A')
[p1,h,STATS] = signrank(res_to_plot(:,1), 0.5, 'tail', 'right')
%calculate effect size
mes(res_to_plot(:,1),0.5,'U3_1')

disp('S1B')
[p2,h,STATS] = signrank(res_to_plot(:,2), 0.5, 'tail', 'left')
%calculate effect size
mes(res_to_plot(:,2),0.5,'U3_1')

disp('S3A')
[p3,h,STATS] = signrank(res_to_plot(:,3), 0.5, 'tail', 'left')
%calculate effect size
mes(res_to_plot(:,3),0.5,'U3_1')

disp('S3B')
[p3,h,STATS] = signrank(res_to_plot(:,4), 0.5, 'tail', 'right')
%calculate effect size
mes(res_to_plot(:,4),0.5,'U3_1')


%% plot results - Figure 1I

cols.k = [0 0 0];
cols.b = [0   15 175];
cols.y = [255 211 0];

color_mat = [cols.k; cols.k; cols.y; cols.k; cols.b; cols.k];

cols.k = [0 0 0];
cols.b = [0 .058 .686];
cols.y = [1  .828 0];
cols.grey = [0.7843 0.7843 0.7843];
cols.dgrey = [0.1922 0.2000 0.2078];

positions = [.25 .75];
pos = [positions-.05; ...
    positions+.05];
   
critical_choices = [res_to_plot(:,1) res_to_plot(:,3)];

figure();
x = 0:1;
plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
hold all;
boxplot(critical_choices, {reshape(['A','B'],2,1) [1;2]}, 'positions',positions, 'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.15,'FactorGap',1, 'Whisker',0); hold all;
scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 200, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',1.5); hold all;
scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,3), 200, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',1.5);
LabelsCS ={'CS^{0}_{.8} vs CS^{0}_{.2}', 'CS^{+}_{.8} vs CS^{+}_{.2}'};
ylim([-.05 1.05]);
xlim([0 1]);
box off
set(findobj(gca,'type','line'),'linew',5)
set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
ybounds = ylim;
set(gca,'YTick',0:0.25:1, 'FontSize',30,'FontName', 'Arial');
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
