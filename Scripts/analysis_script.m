%Analysis script for behavioral and brain-behavioral correlation data (Luettgau, Tempelmann, Kaiser &
%Jocham)

clc; clear all; close all;

sample = 4; %determine experiment to be investigated here: 1 = exp1, 2 = exp2, 3 = exp3, 4 = fmri

exclude_outliers = 1; if exclude_outliers; choice_crit = .50; end %0 = all, 1 = no outliers, 2 = only outliers 

data_path = '...' %set path to exp1, exp 2, exp3 or fmri here

%set path to exp1, exp 2, exp3 or fmri
if sample == 1
   data_path = fullfile([data_path, 'exp1/']);
elseif sample == 2
   data_path = fullfile([data_path, 'exp2/']);
elseif sample == 3
   data_path = fullfile([data_path, 'exp3/']);
elseif sample == 4
    data_path_fmri = fullfile([data_path, 'extracted_parameter_estimates/']);
    data_path = fullfile([data_path, 'fmri/']);
end

cd(data_path)

if sample == 1
   data = dir('candyman*.mat');
elseif sample == 2
   data = dir('candyman2*.mat');
elseif sample == 3
   data = dir('candyman3*.mat');
elseif sample == 4
   data = dir('candyman_fmriLOG_9*.mat');

    %load extracted parameter estimates from fMRI ROIs for brain-behavioral
    %correlations
    
    %Hippocampus
    hippocampus_post_conjunction = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_post_conjunction_anatomical_meants.txt']));
    hippocampus_pre_conjunction = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_pre_conjunction_anatomical_meants.txt']));

    hipp_to_plot = [hippocampus_pre_conjunction -1*(hippocampus_post_conjunction)];
    
    hippocampus_post_2Avs2B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_post_2Avs2B_anatomical_meants.txt']));
    hippocampus_pre_2Avs2B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_pre_2Avs2B_anatomical_meants.txt']));

    hippocampus_post_3Avs3B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_post_3Avs3B_anatomical_meants.txt']));
    hippocampus_pre_3Avs3B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_pre_3Avs3B_anatomical_meants.txt']));
    
    %invert parameter estimates so that they represent associative strength (more positive = more associative strength), easier for interpretation
    inv_hippocampus_post_3Avs3B = -1*(hippocampus_post_3Avs3B); 
    inv_hippocampus_pre_3Avs3B = -1*(hippocampus_pre_3Avs3B);
    
    hippocampus_post_1Avs1B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_post_1Avs1B_anatomical_meants.txt']));
    hippocampus_pre_1Avs1B = load(fullfile([data_path_fmri,'Left_Hippocampus_seclvl_pre_1Avs1B_anatomical_meants.txt']));

    %OFC
    right_ofc_pre_2Avs2B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_pre_2Avs2B_meants.txt']));
    right_ofc_post_2Avs2B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_post_2Avs2B_meants.txt']));

    right_ofc_pre_3Avs3B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_pre_3Avs3B_meants.txt']));
    right_ofc_post_3Avs3B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_post_3Avs3B_meants.txt']));

    right_ofc_pre_1Avs1B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_pre_1Avs1B_meants.txt']));
    right_ofc_post_1Avs1B = load(fullfile([data_path_fmri,'Right_OFC_SVC_seclvl_post_1Avs1B_meants.txt']));

    right_ofc_pre_conjunction = load(fullfile([data_path_fmri,'Right_LOFC_pre_meants.txt']));
    right_ofc_post_conjunction = load(fullfile([data_path_fmri,'Right_LOFC_post_meants.txt']));

    %recode so that positive value is always stronger association and negative value is weaker ass
    ofc_asso_strength = [right_ofc_pre_1Avs1B right_ofc_post_1Avs1B right_ofc_pre_2Avs2B right_ofc_post_2Avs2B  -1*(right_ofc_pre_3Avs3B) -1*(right_ofc_post_3Avs3B)];

    %ofc conjunction
    ofc_to_plot = [right_ofc_pre_conjunction -1*(right_ofc_post_conjunction)];

    %cached value effects
    hippocampus_post_2Avs2B_cachedvalue = load(fullfile([data_path_fmri,'Right_Hippocampus_seclvl_post_cachedvalue_2Avs2B_anatomical_meants.txt']));
    hippocampus_pre_2Avs2B_cachedvalue  = load(fullfile([data_path_fmri,'Right_Hippocampus_seclvl_pre_cachedvalue_2Avs2B_anatomical_meants.txt']));

    hippo_cached_value = [hippocampus_pre_2Avs2B_cachedvalue hippocampus_post_2Avs2B_cachedvalue];
    
    vta_post_2Avs2B_cachedvalue = load(fullfile([data_path_fmri,'Left_VTA_seclvl_post_cachedvalue_2Avs2B_functional_meants.txt']));
    vta_pre_2Avs2B_cachedvalue  = load(fullfile([data_path_fmri,'Left_VTA_seclvl_pre_cachedvalue_2Avs2B_functional_meants.txt']));

    vta_cached_value = [vta_pre_2Avs2B_cachedvalue vta_post_2Avs2B_cachedvalue];
    

    %rois from conjunction contrast thresholded at lenient Z > 2.8
    %(uncorrected)
    lofc_post_3A_3B = load(fullfile([data_path_fmri,'Right_LOFC_seclvl_post_3Avs3B_meants.txt']));
    lofc_post_2A_2B = load(fullfile([data_path_fmri,'Right_LOFC_seclvl_post_2Avs2B_meants.txt']));

    vta_pre_1A_1B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_pre_1Avs1B_meants.txt']));
    vta_post_1A_1B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_post_1Avs1B_meants.txt']));
    
    vta_pre_2A_2B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_pre_2Avs2B_meants.txt']));
    vta_post_2A_2B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_post_2Avs2B_meants.txt']));
    
    vta_pre_3A_3B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_pre_3Avs3B_meants.txt']));
    vta_post_3A_3B = load(fullfile([data_path_fmri,'Left_VTA_seclvl_post_3Avs3B_meants.txt']));
    
    
    accumbens_pre_1A_1B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_pre_1Avs1B_meants.txt']));
    accumbens_post_1A_1B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_post_1Avs1B_meants.txt']));
    
    accumbens_pre_2A_2B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_pre_2Avs2B_meants.txt']));
    accumbens_post_2A_2B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_post_2Avs2B_meants.txt']));
    
    accumbens_pre_3A_3B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_pre_3Avs3B_meants.txt']));
    accumbens_post_3A_3B = load(fullfile([data_path_fmri,'Right_Accumbens_seclvl_post_3Avs3B_meants.txt']));
    

end

format compact 
if exclude_outliers == 1
   exclude = 940; %participant broke up MRI session after first RS run
end


for i = 1:length(data)
    %% load data set
    load(data(i).name)
    
    if sample == 1
       dataset = candyman;
    elseif sample == 2
       dataset = candyman2;
    elseif sample == 3
       dataset = candyman3;
    elseif sample == 4
       dataset = candyman_fmri;
       name = dataset.infoP;
       load(['candyman_fmri_preLOG_',num2str(name),'.mat'])
       load(['candyman_fmri_postLOG_',num2str(name),'.mat'])
       dataset_pre = candyman_fmri_pre;
       dataset_post = candyman_fmri_post;
    end
    
    name = dataset.infoP;
    if name == 999
       name = 162; %due to an error, the wrong subID was logged for this subject
    end
    
    results(i,1) = name;

    if sample == 4 
        if isfield(dataset,'payout')
           results(i,144) = dataset.payout;   
        else
           results(i,144) = NaN;
        end
    end
    
    

    %% Manipulation checks -  button presses for attentional control task during Pavlovian learning
    
    %how many button presses during Pavlovian learning (cover task) 
    if sum(dataset.FOC(:,3) == 0) > length(dataset.FOC(:,3))-(length(dataset.FOC(:,3))*.10)
        results(i,145) = 1; %used for later exlcusion if 1
    end
        results(i,146) = sum(dataset.FOC(:,3) == 0);

    results(:,147) = results(:,1); %repeat subID
    
    %% Manipulation checks - choice-induced revaluation (called devaluation)
    
    %find and delete NaNs (no responses) from choice-induced revaluation trials
    %(called "devaluation")
    %for manipulation check, only trials during choice-induced
    %revaluation trials in which participants saw the learned CS are
    %considered, rest of trials was lures
    if sample ~= 3
        posnan = find(dataset.devaluation(:,2) == 1 | dataset.devaluation(:,3) == 1 | ...
                      dataset.devaluation(:,2) == 2 | dataset.devaluation(:,3) == 2); 
    else
        posnan = find(dataset.devaluation(:,2) == 1 | dataset.devaluation(:,2) == 2 | ...
                      dataset.devaluation(:,2) == 3); 
    end
    devaluation = dataset.devaluation;
    devaluation = devaluation(posnan,:);
    devalpos = find(isnan(devaluation(:,12)));
    %delete NaNs and other trials of no interest
    devaluation(devalpos,:) = [];

    %define choice R/L, where the CS associated with the higher valued US was
    %chosen
    highvalchoiceL = length(find(devaluation(:,2) > devaluation(:,3) & devaluation(:,12) == 1));
    highvalchoiceR = length(find(devaluation(:,3) > devaluation(:,2) & devaluation(:,12) == 2));

    %compute percent choices of higher valued CS during choice-induced
    %revaluation
    percentchoicehighval = (highvalchoiceL + highvalchoiceR)/length(devaluation);
    results(i,2) = percentchoicehighval;
    
    
    if sample == 3
       %calculate choice prob. between CS-A and CS0A
       %terminology is CS1A = CS-A, CS1B = CS-B, CS2A = CS0A, CS2B = CS0B, CS3A = CS+A, CS3B = CS+B
       CS2AvsCS1AchoiceL = length(find(devaluation(:,2) == 2 & devaluation(:,3) == 1 & devaluation(:,12) == 1));
       CS2AvsCS1AchoiceR = length(find(devaluation(:,2) == 1 & devaluation(:,3) == 2 & devaluation(:,12) == 2));
       percentchoiceCS2AvsCS1A = (CS2AvsCS1AchoiceL+CS2AvsCS1AchoiceR)/(length(find(devaluation(:,2) == 1 & devaluation(:,3) == 2)) + length(find(devaluation(:,2) == 2 & devaluation(:,3) == 1)));
       results(i,143) = percentchoiceCS2AvsCS1A;
        
       %calculate choice prob. between CS0A and CS+A 
       CS3AvsCS2AchoiceL = length(find(devaluation(:,2) == 3 & devaluation(:,3) == 2 & devaluation(:,12) == 1));
       CS3AvsCS2AchoiceR = length(find(devaluation(:,2) == 2 & devaluation(:,3) == 3 & devaluation(:,12) == 2));
       percentchoiceCS3AvsCS2A = (CS3AvsCS2AchoiceL +CS3AvsCS2AchoiceR)/(length(find(devaluation(:,2) == 3 & devaluation(:,3) == 2)) + length(find(devaluation(:,2) == 2 & devaluation(:,3) == 3)));
       results(i,144) = percentchoiceCS3AvsCS2A;

    end
               
    %% Decision Choice Probe (called forcedchoicekanjis)

    %delete NaNs from forced choice responses (kanjis)
    posnankanjis = find(isnan(dataset.forcedchoicekanjis(:,12)));
    forcedchoicekanjis = dataset.forcedchoicekanjis;
    forcedchoicekanjis(posnankanjis,:) = [];

    %define choice R/L, where the CS associated with the higher valued US was
    %chosen - kanjis
    highvalkanjichoiceL = length(find(forcedchoicekanjis(:,7) > forcedchoicekanjis(:,11) & forcedchoicekanjis(:,12) == 1));
    highvalkanjichoiceR = length(find(forcedchoicekanjis(:,11) > forcedchoicekanjis(:,7) & forcedchoicekanjis(:,12) == 2));

    %percent choices of higher valued CS
    percentchoicehighvalkanji = (highvalkanjichoiceL + highvalkanjichoiceR)/(length(find(forcedchoicekanjis(:,11) > forcedchoicekanjis(:,7))) + length(find(forcedchoicekanjis(:,7) > forcedchoicekanjis(:,11))));
    results(i,3) = percentchoicehighvalkanji;


    %% Decision Choice Probe - loop over all CS for all binary combinations of choices and compute choice probability
    for m = 1:6 %order is 1 = CS-A, 2 = CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
        [posL] = find(forcedchoicekanjis(:,2) == m);
        [posR] = find(forcedchoicekanjis(:,3) == m);

        CS(m) = (sum(forcedchoicekanjis(posL,12) == 1) + sum(forcedchoicekanjis(posR,12) == 2)) / (length(posL) + length(posR));
    end
    
    %order is CS-A, CS-B, CS0A, CS0B, CS+A, CS+B
    results(i,4:9) = CS;
    
    %% Decision Choice Probe - loop over all single value categroy specific comparisons (supplementary figure 2)
    %m = chosen, n = unchosen
    %order is 1 = CS-A, 2 = CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
    for m = 1:6
        for n = 1:6
            if n ~= m 
                CSvCS = sum(forcedchoicekanjis(:,2) == m & forcedchoicekanjis(:,3) == n & forcedchoicekanjis(:,12) == 1 | forcedchoicekanjis(:,2) == n & forcedchoicekanjis(:,3) == m & forcedchoicekanjis(:,12) == 2);
                percentchoiceCS = CSvCS/sum(forcedchoicekanjis(:,2) == m & forcedchoicekanjis(:,3) == n | forcedchoicekanjis(:,2) == n & forcedchoicekanjis(:,3) == m);
                pairwise_comp(n,m) = percentchoiceCS;
            end
        end
    end

    %step through matrix and extract pairwise choice probabilities
    elems = [2 1; 3 1; 4 1; 5 1; 6 1; 3 2; 4 2; 5 2; 6 2; 4 3; 5 3; 6 3; 5 4; 6 4; 6 5];
    count = 1;
    for k = 10:2:38
        results(i,k) = pairwise_comp(elems(count,1), elems(count,2));
        results(i,k+1) = 1-results(i,k);
        count = count+1;
    end


%% Difference score pre- to post-exp rating of CS

if sample == 4
    for j = 1:6
        prepostkanjis = find(dataset_pre.ratingkanjis(5,j) == dataset_post.postratingkanjis(5,:));
        value = find(dataset_pre.ratingkanjis(5,j) == dataset_pre.FOC(:,5));
        valuekanji = dataset_pre.FOC(value(1),8) ;

        diffkanji(1,j) = dataset_pre.ratingkanjis(2,j);
        diffkanji(2,j) = dataset_post.postratingkanjis(2,prepostkanjis);
        diffkanji(3,j) = dataset_post.postratingkanjis(2,prepostkanjis) - dataset_pre.ratingkanjis(2,j);
        diffkanji(4,j) = valuekanji;
        diffkanji(5,j) = dataset_pre.FOC(value(1),5);
    end
    diffkanji = diffkanji';
    diffkanji = sortrows(diffkanji,4);
    diffkanji = diffkanji';
    results(i,40:45) = diffkanji(1,:);

    results(i,46:51) = diffkanji(2,:);

    results(i,52:57) = diffkanji(3,:);

    results(i,58:63) = diffkanji(4,:);
    
    j = [];
    

    %% Difference score pre- to post-exp rating of US

    for j = 1:3

        val =  unique(dataset_pre.FOC(:,7));
        prefood = find(dataset_pre.ratingfoodcues(5,:) == val(j));
        postfood = find(dataset_post.postratingfoodcues(5,:) == val(j));
        valfood = dataset_pre.FOC(val(j),8);

        difffood(1,j) = val(j);
        difffood(2,j) = dataset_pre.ratingfoodcues(2,prefood);
        difffood(3,j) = dataset_post.postratingfoodcues(2,postfood);
        difffood(4,j) = dataset_post.postratingfoodcues(2,postfood) - dataset_pre.ratingfoodcues(2,prefood);

    end
    difffood = difffood';
    difffood = sortrows(difffood,2);
    difffood = difffood';
    results(i,64:66) = difffood(1,:);

    results(i,67:69) = difffood(2,:);

    results(i,70:72) = difffood(3,:);

    results(i,73:75) = difffood(4,:);

    
else
    
    
        for j = 1:6
            prepostkanjis = find(dataset.ratingkanjis(5,j) == dataset.postratingkanjis(5,:));
            value = find(dataset.ratingkanjis(5,j) == dataset.FOC(:,5));
            valuekanji = dataset.FOC(value(1),8) ;

            diffkanji(1,j) = dataset.ratingkanjis(2,j);
            diffkanji(2,j) = dataset.postratingkanjis(2,prepostkanjis);
            diffkanji(3,j) = dataset.postratingkanjis(2,prepostkanjis) - dataset.ratingkanjis(2,j);
            diffkanji(4,j) = valuekanji;
            diffkanji(5,j) = dataset.FOC(value(1),5);
        end
        diffkanji = diffkanji';
        diffkanji = sortrows(diffkanji,4);
        diffkanji = diffkanji';
        results(i,40:45) = diffkanji(1,:);

        results(i,46:51) = diffkanji(2,:);

        results(i,52:57) = diffkanji(3,:);

        results(i,58:63) = diffkanji(4,:);

        j = [];

        %% Difference score pre- to post-exp rating of US
    
        for j = 1:3

            val =  unique(dataset.FOC(:,7));
            prefood = find(dataset.ratingfoodcues(5,:) == val(j));
         if isfield(dataset,'postratingfoodcues')   
            postfood = find(dataset.postratingfoodcues(5,:) == val(j));
         end
            valfood = dataset.FOC(val(j),8);

            difffood(1,j) = val(j);
            difffood(2,j) = dataset.ratingfoodcues(2,prefood);
         if isfield(dataset,'postratingfoodcues')   
            difffood(3,j) = dataset.postratingfoodcues(2,postfood);
            difffood(4,j) = dataset.postratingfoodcues(2,postfood) - dataset.ratingfoodcues(2,prefood);
         end

        end
        difffood = difffood';
        difffood = sortrows(difffood,2);
        difffood = difffood';
        results(i,64:66) = difffood(1,:);

        results(i,67:69) = difffood(2,:);
      
    if isfield(dataset,'postratingfoodcues')   
        results(i,70:72) = difffood(3,:);
        results(i,73:75) = difffood(4,:);
    else
        results(i,70:75) = NaN;
    end

end




%% responses during the attentional control task in fMRI
if sample == 4
    
    prereptrials = find(dataset.prerepsup(:,11) == 1);
    repsup_pre = dataset.prerepsup(prereptrials,:);
    
    postreptrials = find(dataset.postrepsup(:,11) == 1);
    repsup_post = dataset.postrepsup(postreptrials,:);
    
    results(i,142) = sum(repsup_pre(:,18)) + sum(repsup_post(:,18)) ;
    results(i,143) = results(i,142)/144;

    %correct answers during rep sup pre

    repsup_pre = repsup_pre((~isnan(repsup_pre(:,15)) == 1),:);

    correct_1_1_trials_pre = length(find(repsup_pre(:,2) == 1 & repsup_pre(:,3) == 1 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 1 & repsup_pre(:,3) == 1));
    results(i,84) = correct_1_1_trials_pre;

    correct_2_1_trials_pre = length(find(repsup_pre(:,2) == 2 & repsup_pre(:,3) == 1 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 2 & repsup_pre(:,3) == 1));
    results(i,87) = correct_2_1_trials_pre;

    correct_3_2_trials_pre = length(find(repsup_pre(:,2) == 3 & repsup_pre(:,3) == 2 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 3 & repsup_pre(:,3) == 2));
    results(i,90) = correct_3_2_trials_pre;

    correct_4_2_trials_pre = length(find(repsup_pre(:,2) == 4 & repsup_pre(:,3) == 2 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 4 & repsup_pre(:,3) == 2));
    results(i,93) = correct_4_2_trials_pre;

    correct_5_3_trials_pre = length(find(repsup_pre(:,2) == 5 & repsup_pre(:,3) == 3 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 5 & repsup_pre(:,3) == 3));
    results(i,96) = correct_5_3_trials_pre;

    correct_6_3_trials_pre = length(find(repsup_pre(:,2) == 6 & repsup_pre(:,3) == 3 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 6 & repsup_pre(:,3) == 3));
    results(i,99) = correct_6_3_trials_pre;

    correct_1_trials_pre = length(find(repsup_pre(:,2) == 1 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 1));
    results(i,106) = correct_1_trials_pre;

    correct_2_trials_pre = length(find(repsup_pre(:,2) == 2 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 2));
    results(i,107) = correct_2_trials_pre;

    correct_3_trials_pre = length(find(repsup_pre(:,2) == 3 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 3));
    results(i,108) = correct_3_trials_pre;

    correct_4_trials_pre = length(find(repsup_pre(:,2) == 4 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 4));
    results(i,109) = correct_4_trials_pre;

    correct_5_trials_pre = length(find(repsup_pre(:,2) == 5 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 5));
    results(i,110) = correct_5_trials_pre;

    correct_6_trials_pre = length(find(repsup_pre(:,2) == 6 & repsup_pre(:,18) == 1))/length(find(repsup_pre(:,2) == 6));
    results(i,111) = correct_6_trials_pre;

    %correct answers during rep sup post
    repsup_post = repsup_post((~isnan(repsup_post(:,15)) == 1),:);

    correct_1_1_trials_post = length(find(repsup_post(:,2) == 1 & repsup_post(:,3) == 1 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 1 & repsup_post(:,3) == 1));
    results(i,85) = correct_1_1_trials_post;
    results(i,86) = (results(i,84) + results(i,85))/2;

    correct_2_1_trials_post = length(find(repsup_post(:,2) == 2 & repsup_post(:,3) == 1 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 2 & repsup_post(:,3) == 1));
    results(i,88) = correct_2_1_trials_post;
    results(i,89) = (results(i,87) + results(i,88))/2;

    correct_3_2_trials_post = length(find(repsup_post(:,2) == 3 & repsup_post(:,3) == 2 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 3 & repsup_post(:,3) == 2));
    results(i,91) = correct_3_2_trials_post;
    results(i,92) = (results(i,90) + results(i,91))/2;

    correct_4_2_trials_post = length(find(repsup_post(:,2) == 4 & repsup_post(:,3) == 2 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 4 & repsup_post(:,3) == 2));
    results(i,94) = correct_4_2_trials_post;
    results(i,95) = (results(i,93) + results(i,94))/2;

    correct_5_3_trials_post = length(find(repsup_post(:,2) == 5 & repsup_post(:,3) == 3 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 5 & repsup_post(:,3) == 3));
    results(i,97) = correct_5_3_trials_post;
    results(i,98) = (results(i,96) + results(i,97))/2;

    correct_6_3_trials_post = length(find(repsup_post(:,2) == 6 & repsup_post(:,3) == 3 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 6 & repsup_post(:,3) == 3));
    results(i,100) = correct_6_3_trials_post;
    results(i,101) = (results(i,99)+results(i,100))/2;

    results(i,102) = sum(results(i,86:3:101))/6;


    correct_1_trials_post = length(find(repsup_post(:,2) == 1 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 1));
    results(i,112) = correct_1_trials_post;

    correct_2_trials_post = length(find(repsup_post(:,2) == 2 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 2));
    results(i,113) = correct_2_trials_post;

    correct_3_trials_post = length(find(repsup_post(:,2) == 3 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 3));
    results(i,114) = correct_3_trials_post;

    correct_4_trials_post = length(find(repsup_post(:,2) == 4 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 4));
    results(i,115) = correct_4_trials_post;

    correct_5_trials_post = length(find(repsup_post(:,2) == 5 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 5));
    results(i,116) = correct_5_trials_post;

    correct_6_trials_post = length(find(repsup_post(:,2) == 6 & repsup_post(:,18) == 1))/length(find(repsup_post(:,2) == 6));
    results(i,117) = correct_6_trials_post;

    results(i,118) = sum(results(i,106:111))/6;%average correct pre
    results(i,119) = sum(results(i,112:117))/6;%average correct post

    results(i,120) = (results(i,118)+results(i,119))/2;%average correct pre and post
end


%automated exclusion of outliers based on choice-induced revaluation
%behavior

if exclude_outliers == 1
    if results(i,2) < choice_crit
       results(i,2:end) = NaN;
       range_kanjis(i,:) = NaN;
    end
elseif exclude_outliers == 2
    if results(i,2) > choice_crit
       results(i,2:end) = NaN;
    end
end


if exclude_outliers == 1
   if sum(results(i,1) == exclude) > 0
      results(i,2:end) = NaN;
   end
end

end

%automated exclusion based on non-responses during Pavlovian learning phase cover task
pos_rt = find(results(:,145) == 0);
results = results(pos_rt,:);


if exclude_outliers == 1
    results = results(~isnan(results(:,2)),:);
elseif exclude_outliers == 2
    results = results(~isnan(results(:,2)),:);   
end

%% statistics

%% descriptive statistics - US subjective value (pre-rating)
pre_food = results(:,67:69);
nanmean(pre_food)
nanstd(pre_food)

%attentional control task performance (fMRI)
if sample == 4
    disp('Percent correct answers in Repetition Suppression overall')
    mean(results(:,120))
    std(results(:,120))
    [h,p,ci,stats] = ttest(results(:,120), 0.5)


    disp('Percent correct answers in Repetition Suppression pre')
    mean(results(:,118))
    std(results(:,118))
    [h,p,ci,stats] = ttest(results(:,118), 0.5)

    disp('Percent correct answers in Repetition Suppression post')
    mean(results(:,119))
    std(results(:,119))
    [h,p,ci,stats] = ttest(results(:,119), 0.5)


    [h,p,ci,stats] = ttest(results(:,118), results(:,119))

    mean(results(:,144))
    std(results(:,144))
end
%% CS subjective value (pre-rating)   
pre_kanji = results(:,40:45);
nanmean(pre_kanji)
nanstd(pre_kanji)

% repeated measures ANOVA on CS subjective value (pre-rating) 
data = pre_kanji;
varNames = {'S1','S2','S3','S4','S5','S6'};

t = array2table(data,'VariableNames',varNames);

factorNames = {'stimulus'};
within = table({'one'; 'two';'three';'four';'five';'six'},'VariableNames',factorNames);


rm = fitrm(t,'S1-S6~1','WithinDesign',within); 


[ranovatbl] = ranova(rm, 'WithinModel', 'stimulus');

disp('rmANOVA on kanjis ratings')
ranovatbl

%% overall choice probabilities for CS during choice probe phase
%terminology is CS1A = CS-A, CS1B = CS-B, CS2A = CS0A, CS2B = CS0B, CS3A = CS+A, CS3B = CS+B

CS1A = results(:,4);
CS1B = results(:,5);
CS2A = results(:,6);
CS2B = results(:,7);
CS3A = results(:,8);
CS3B = results(:,9);

results_kanjis = results(:,10:2:38);

results_kanjis_to_plot = results_kanjis;

res_to_plot = [CS1A CS1B CS2A CS2B CS3A CS3B];

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
[P,H,STATS] = signrank(CS1A, CS1B)

%CS2A = CS0A vs. CS2B = CS0B
disp('overall choice probability CS0A vs. CS0B')
[P,H,STATS] = signrank(CS2A, CS2B)


%CS3A = CS+A vs. CS3B = CS+B
disp('overall choice probability CS+A vs. CS+B')
[P,H,STATS] = signrank(CS3A, CS3B)



%% plot behavioral results
%used for Figure 1E-H

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

if sample == 1
    %exp1
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0); hold all; 
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5); hold all;
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',2.5);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',2.5);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
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


elseif sample == 2
    %exp2
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0);
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',2.5);
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',2.5);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
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


elseif sample == 3
    %exp3
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0);
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',2.5);
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.b,'LineWidth',2.5);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',2.5);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.k, 'MarkerEdgeColor', cols.k,'LineWidth',2.5);
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


elseif sample == 4
    %fmri
    figure();
    x = 0:9;
    plot(x,ones(1,length(x))*0.5, 'k--', 'linewidth', 4);
    hold all; 
    boxplot(res_to_plot, {reshape(repmat('A':'C',2,1),6,1) repmat((1:2)',3,1)},'boxstyle', 'outline', 'colors', cols.k, 'symbol','','Widths',0.6,'FactorGap',20, 'Whisker',0);
    scatter(linspace(pos(1,1), pos(2,1),length(res_to_plot)), res_to_plot(:,1), 150, 'o', 'MarkerFaceColor', cols.dgrey, 'MarkerEdgeColor', cols.dgrey,'LineWidth',2.5); 
    scatter(linspace(pos(1,2), pos(2,2), length(res_to_plot)), res_to_plot(:,2), 150, 'o', 'MarkerFaceColor', cols.dgrey, 'MarkerEdgeColor', cols.dgrey,'LineWidth',2.5);
    scatter(linspace(pos(1,3), pos(2,3), length(res_to_plot)), res_to_plot(:,3), 150, 'o', 'MarkerFaceColor', cols.y, 'MarkerEdgeColor', cols.y,'LineWidth',2.5);
    scatter(linspace(pos(1,4), pos(2,4), length(res_to_plot)), res_to_plot(:,4), 150, 'o', 'MarkerFaceColor', cols.dgrey, 'MarkerEdgeColor', cols.dgrey,'LineWidth',2.5);
    scatter(linspace(pos(1,5), pos(2,5), length(res_to_plot)), res_to_plot(:,5), 150, 'o', 'MarkerFaceColor', cols.b, 'MarkerEdgeColor', cols.b,'LineWidth',2.5);
    scatter(linspace(pos(1,6), pos(2,6), length(res_to_plot)), res_to_plot(:,6), 150, 'o', 'MarkerFaceColor', cols.dgrey, 'MarkerEdgeColor', cols.dgrey,'LineWidth',2.5);
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


end

if sample == 4
    %Tests on extracted fMRI parameter estimates
    
    % Hippocampus
    hipp_asso_strength = [hippocampus_pre_1Avs1B hippocampus_post_1Avs1B hippocampus_pre_2Avs2B hippocampus_post_2Avs2B -1*(hippocampus_pre_3Avs3B) -1*(hippocampus_post_3Avs3B)];
    %recode so that positive value is always stronger association and negative
    %value is weaker association
    
    % repeated measures ANOVA
    data = hipp_asso_strength(:,3:end);
    varNames = {'S2A2Bpre', 'S2A2Bpost', 'S3A3Bpre', 'S3A3Bpost'};

    t = array2table(data,'VariableNames',varNames);

    factorNames = {'valence','time'};
    within = table({'med';'med';'hi';'hi'},{'pre';'post';'pre';'post'},'VariableNames',factorNames);


    rm = fitrm(t,'S2A2Bpre-S3A3Bpost~1','WithinDesign',within); 


    [ranovatbl] = ranova(rm, 'WithinModel', 'valence*time');

    [P,H,STATS] = signrank(hipp_asso_strength(:,1), 0)
    [P,H,STATS] = signrank(hipp_asso_strength(:,2), 0)
    [P,H,STATS] = signrank(hipp_asso_strength(:,3), 0)
    [P,H,STATS] = signrank(hipp_asso_strength(:,4), 0)
    [P,H,STATS] = signrank(hipp_asso_strength(:,5), 0)
    [P,H,STATS] = signrank(hipp_asso_strength(:,6), 0)
    
    [P,H,STATS] = signrank(hipp_asso_strength(:,1), hipp_asso_strength(:,2))
    [P,H,STATS] = signrank(hipp_asso_strength(:,3), hipp_asso_strength(:,4))
    [P,H,STATS] = signrank(hipp_asso_strength(:,5), hipp_asso_strength(:,6))
    
    %OFC

    % repeated measures ANOVA
    data = ofc_asso_strength(:,3:end);
    varNames = {'S2A2Bpre', 'S2A2Bpost', 'S3A3Bpre', 'S3A3Bpost'};

    t = array2table(data,'VariableNames',varNames);

    factorNames = {'valence','time'};
    within = table({'med';'med';'hi';'hi'},{'pre';'post';'pre';'post'},'VariableNames',factorNames);


    rm = fitrm(t,'S2A2Bpre-S3A3Bpost~1','WithinDesign',within); 


    [ranovatbl] = ranova(rm, 'WithinModel', 'valence*time');


    [P,H,STATS] = signrank(ofc_asso_strength(:,1), 0)
    [P,H,STATS] = signrank(ofc_asso_strength(:,2), 0)
    [P,H,STATS] = signrank(ofc_asso_strength(:,3), 0)
    [P,H,STATS] = signrank(ofc_asso_strength(:,4), 0)
    [P,H,STATS] = signrank(ofc_asso_strength(:,5), 0)
    [P,H,STATS] = signrank(ofc_asso_strength(:,6), 0)

    [P,H,STATS] = signrank(ofc_asso_strength(:,1), ofc_asso_strength(:,2))
    [P,H,STATS] = signrank(ofc_asso_strength(:,3), ofc_asso_strength(:,4))
    [P,H,STATS] = signrank(ofc_asso_strength(:,5), ofc_asso_strength(:,6))

    
    %% used this plot for illustration in paper - Figure 2D&E
    
    cols.k = [0 0 0];
    cols.b = [0 .058 .686];
    cols.y = [1  .828 0];
    cols.grey = [0.7843 0.7843 0.7843];
    cols.dgrey = [0.1922 0.2000 0.2078];


    b = figure();
    bar(nanmean(hipp_to_plot), 'k', 'BarWidth', 0.6); hold all;
    errorbar(1:2,nanmean(hipp_to_plot),std(hipp_to_plot)./sqrt(size(hipp_to_plot,1)), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 2.5);
    Labels ={'PRE', 'POST'};
    ylim([-50 100]); 
    xlim([0.5 2.5]); 
    set(gca,'XTick',1:2,'XTickLabel', Labels);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):50:ybounds(2), 'FontSize',15,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'XTickLabel', Labels, 'FontSize',15,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    ticklabels_new = cell(size(Labels));
    for i = 1:length(Labels)
        ticklabels_new{i} = ['\color{black} ' Labels{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new,'FontSize',15,'FontName', 'Arial');
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',15,'FontName', 'Arial');

    positions = [1 2 4 5 7 8];   
    pos = [positions-1; ...
           positions+1];

    mean_hipp_asso_strength = [mean(hippocampus_pre_1Avs1B);
                               mean(hippocampus_post_1Avs1B);
                               NaN;
                               mean(hippocampus_pre_2Avs2B);
                               mean(hippocampus_post_2Avs2B);
                               NaN;
                               mean(-1*(hippocampus_pre_3Avs3B));
                               mean(-1*(hippocampus_post_3Avs3B))];
    size_vec = [1 2 4 5 7 8];
    color_vec = [cols.k; cols.k; cols.k; cols.y; cols.y; cols.k; cols.b; cols.b];

    a = figure();  
    bar(1:8,mean_hipp_asso_strength, 'k', 'BarWidth', 0.6); hold all;
    errorbar(size_vec,mean_hipp_asso_strength(size_vec,:), std(hipp_asso_strength)./sqrt(size(hipp_asso_strength,1)), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 2.5);
    LabelsCS ={'PRE', 'POST','PRE', 'POST', 'PRE', 'POST'};
    ylim([-80 80]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):40:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',size_vec)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new,'FontSize',30,'FontName', 'Arial');
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',30,'FontName', 'Arial');

%     inset(a, b, .20);
%     set(gcf,'color','w');

    
    c = figure();
    bar(nanmean(ofc_to_plot), 'k', 'BarWidth', 0.6); hold all;
    errorbar(1:2,nanmean(ofc_to_plot),std(ofc_to_plot)./sqrt(size(ofc_to_plot,1)), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 2.5);
    Labels ={'PRE', 'POST'};
    ylim([-50 100]); 
    xlim([0.5 2.5]); 
    set(gca,'XTick',1:2,'XTickLabel', Labels);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):50:ybounds(2), 'FontSize',15,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'XTickLabel', Labels, 'FontSize',15,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    ticklabels_new = cell(size(Labels));
    for i = 1:length(Labels)
        ticklabels_new{i} = ['\color{black} ' Labels{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new,'FontSize',15,'FontName', 'Arial');
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',15,'FontName', 'Arial');



    mean_ofc_asso_strength = [mean(right_ofc_pre_1Avs1B);
                               mean(right_ofc_post_1Avs1B);
                               NaN;
                               mean(right_ofc_pre_2Avs2B);
                               mean(right_ofc_post_2Avs2B);
                               NaN;
                               mean(-1*(right_ofc_pre_3Avs3B));
                               mean(-1*(right_ofc_post_3Avs3B))];

    d = figure();  
    bar(1:8,mean_ofc_asso_strength, 'k', 'BarWidth', 0.6); hold all;
    errorbar(size_vec,mean_ofc_asso_strength(size_vec,:), std(ofc_asso_strength)./sqrt(size(ofc_asso_strength,1)), 'linestyle', 'none', 'color', 'k', 'CapSize', 0, 'LineWidth', 2.5);
    LabelsCS ={'PRE', 'POST','PRE', 'POST', 'PRE', 'POST'};
    ylim([-80 80]); 
    xlim([0 9]);
    box off
    set(findobj(gca,'type','line'),'linew',5)
    set(gca,'TickLength',[0.01, 0.001],'linewidth',2.5)
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):40:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'xtick',size_vec)
    set(gca,'XTickLabel', LabelsCS, 'FontSize',30,'FontName', 'Arial');
    set(gca,'XTickLabelRotation', 45);
    set(gcf,'color','w');
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    ticklabels_new = cell(size(LabelsCS));
    for i = 1:length(LabelsCS)
        ticklabels_new{i} = ['\color{black} ' LabelsCS{i}];
    end
    % set the tick labels
    set(gca, 'XTickLabel', ticklabels_new,'FontSize',30,'FontName', 'Arial');
    % prepend a color for each tick label
    LabelsY = get(gca,'YTickLabel');
    ticklabels_ynew = cell(size(LabelsY));
    for i = 1:length(LabelsY)
        ticklabels_ynew{i} = ['\color{black} ' LabelsY{i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_ynew,'FontSize',30,'FontName', 'Arial');


%     inset(d, c, .20);
%     set(gcf,'color','w');




    
    %% used this plot for illustration in paper - Figure 3B

    %calculation of overall choice prob. differences
    diff_3A_3B = (res_to_plot(:,5) - res_to_plot(:,6));
    diff_2A_2B = (res_to_plot(:,3) - res_to_plot(:,4));
    
    
    %hippocampus cluster asso strength from post seclvl 3Avs3B - correlation with choice
    %behavior


    [rho,p] = corr(diff_3A_3B,inv_hippocampus_post_3Avs3B, 'type', 'Spearman')


    cols.k = [0 0 0];
    cols.b = [0 .058 .686];
    cols.y = [1  .828 0];
    cols.grey = [0.7843 0.7843 0.7843];

    figure;
    scatter(diff_3A_3B, inv_hippocampus_post_3Avs3B, 200, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    xlim([-.4 1]); 
    xticks(-.4:.4:1);
    ylim([-100 200]);
    ybounds = ylim;
    box off
    regline = lsline;  
    regline.LineWidth = 2.5;  
    regline.Color = cols.k;
    regline.XData = [-.3 .86];
    regline.LineStyle = '--';
    set(gca,'FontSize',30,'FontName', 'Arial');   
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    set(gca,'YTick',ybounds(1):100:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    set(gcf,'color','w');



    %% use this plot for illustration in paper - Figure S2A
    %calculation of within-category choice prob. differences
    results_pairs = results_kanjis_to_plot(:,[1,10,15]);

    %hippocampus cluster asso strength from post seclvl 3Avs3B - correlation with choice
    %behavior

    [rho,p] = corr(results_pairs(:,3),inv_hippocampus_post_3Avs3B, 'type', 'Spearman') 

    cols.k = [0 0 0];
    cols.b = [0 .058 .686];
    cols.y = [1  .828 0];
    cols.grey = [0.7843 0.7843 0.7843];

    figure;
    scatter(results_pairs(:,3), inv_hippocampus_post_3Avs3B, 200, 'filled', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0]);
    xlim([0 1]); 
    xticks(0:.25:1);
    ylim([-100 200]);
    ybounds = ylim;
    box off
    regline = lsline;  
    regline.LineWidth = 2.5;  
    regline.Color = cols.k;
    regline.XData = [.1 1];
    regline.LineStyle = '--';
    set(gca,'FontSize',30,'FontName', 'Arial');   
    set(gca,'TickLength',[0.01, 0.001],'linewidth',1.5)
    set(gca,'YTick',ybounds(1):100:ybounds(2), 'FontSize',30,'FontName', 'Arial');
    set(gca,'TickDir','out')
    set(gca,'ycolor',cols.k)
    set(gca,'xcolor',cols.k)
    set(gcf,'color','w'); %hold on


    %% sec lvl post ROI (thresh at z=2.8)   
    
    %2A-2B - right LOFC
    [rho,p] = corr(results_pairs(:,2),lofc_post_2A_2B, 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,lofc_post_2A_2B, 'type', 'Spearman')

    %3A-3B - right LOFC
    [rho,p] = corr(diff_3A_3B,lofc_post_3A_3B, 'type', 'Spearman')
    [rho,p] = corr(results_pairs(:,3),lofc_post_3A_3B, 'type', 'Spearman')

    %NAcc
    nacc_asso_strength = [accumbens_pre_1A_1B accumbens_post_1A_1B accumbens_pre_2A_2B accumbens_post_2A_2B -1*(accumbens_pre_3A_3B) -1*(accumbens_post_3A_3B)];
    
    % repeated measures ANOVA
    data = nacc_asso_strength;
    varNames = {'S1A1Bpre', 'S1A1Bpost','S2A2Bpre', 'S2A2Bpost', 'S3A3Bpre', 'S3A3Bpost'};

    t = array2table(data,'VariableNames',varNames);

    factorNames = {'valence','time'};
    within = table({'low';'low';'med';'med';'hi';'hi'},{'pre';'post';'pre';'post';'pre';'post'},'VariableNames',factorNames);


    rm = fitrm(t,'S1A1Bpre-S3A3Bpost~1','WithinDesign',within); 


    [ranovatbl] = ranova(rm, 'WithinModel', 'valence*time');
    
    [P,H,STATS] = signrank(nacc_asso_strength(:,1), 0)
    [P,H,STATS] = signrank(nacc_asso_strength(:,2), 0)
    [P,H,STATS] = signrank(nacc_asso_strength(:,3), 0)
    [P,H,STATS] = signrank(nacc_asso_strength(:,4), 0)
    [P,H,STATS] = signrank(nacc_asso_strength(:,5), 0)
    [P,H,STATS] = signrank(nacc_asso_strength(:,6), 0)

    [P,H,STATS] = signrank(nacc_asso_strength(:,1), nacc_asso_strength(:,2))
    [P,H,STATS] = signrank(nacc_asso_strength(:,3), nacc_asso_strength(:,4))
    [P,H,STATS] = signrank(nacc_asso_strength(:,5), nacc_asso_strength(:,6))
    
    %2A-2B - right Accumbens
    [rho,p] = corr(results_pairs(:,2),accumbens_post_2A_2B, 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,accumbens_post_2A_2B, 'type', 'Spearman')
    
    %3A-3B - right Accumbens
    [rho,p] = corr(diff_3A_3B,accumbens_post_3A_3B, 'type', 'Spearman')
    [rho,p] = corr(results_pairs(:,3),accumbens_post_3A_3B, 'type', 'Spearman')
    

    % repeated measures ANOVA
    vta_asso_strength = [vta_pre_1A_1B vta_post_1A_1B vta_pre_2A_2B vta_post_2A_2B -1*(vta_pre_3A_3B) -1*(vta_post_3A_3B)];
    
    data = vta_asso_strength;
    varNames = {'S1A1Bpre', 'S1A1Bpost','S2A2Bpre', 'S2A2Bpost', 'S3A3Bpre', 'S3A3Bpost'};

    t = array2table(data,'VariableNames',varNames);

    factorNames = {'valence','time'};
    within = table({'low';'low';'med';'med';'hi';'hi'},{'pre';'post';'pre';'post';'pre';'post'},'VariableNames',factorNames);


    rm = fitrm(t,'S1A1Bpre-S3A3Bpost~1','WithinDesign',within); 

    [ranovatbl] = ranova(rm, 'WithinModel', 'valence*time');

    [P,H,STATS] = signrank(vta_asso_strength(:,1), 0)
    [P,H,STATS] = signrank(vta_asso_strength(:,2), 0)
    [P,H,STATS] = signrank(vta_asso_strength(:,3), 0)
    [P,H,STATS] = signrank(vta_asso_strength(:,4), 0)
    [P,H,STATS] = signrank(vta_asso_strength(:,5), 0)
    [P,H,STATS] = signrank(vta_asso_strength(:,6), 0)

    [P,H,STATS] = signrank(vta_asso_strength(:,1), vta_asso_strength(:,2))
    [P,H,STATS] = signrank(vta_asso_strength(:,3), vta_asso_strength(:,4))
    [P,H,STATS] = signrank(vta_asso_strength(:,5), vta_asso_strength(:,6))

    %2A-2B - left VTA
    [rho,p] = corr(results_pairs(:,2),vta_post_2A_2B, 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,vta_post_2A_2B, 'type', 'Spearman')

    %3A-3B - left VTA
    [rho,p] = corr(diff_3A_3B,vta_post_3A_3B, 'type', 'Spearman')
    [rho,p] = corr(results_pairs(:,3),vta_post_3A_3B, 'type', 'Spearman')



    %% cached value control analysis - similarity with US- 

    [P,H,STATS] = signrank(hippo_cached_value(:,1), 0)
    [P,H,STATS] = signrank(hippo_cached_value(:,2), 0)

    [P,H,STATS] = signrank(hippo_cached_value(:,1), hippo_cached_value(:,2))
        
    %correlation of cached-value effect (R Hippocampus) and associative effect (L Hippocampus)
    %post - 2A-2B
    [rho,p] = corr(hippocampus_post_2Avs2B,hippocampus_post_2Avs2B_cachedvalue, 'type', 'Spearman')

    %post - conjunction
    [rho,p] = corr(hippocampus_post_conjunction,hippocampus_post_2Avs2B_cachedvalue, 'type', 'Spearman')

    %pre
    [rho,p] = corr(hippocampus_pre_2Avs2B,hippocampus_pre_2Avs2B_cachedvalue, 'type', 'Spearman')
    
    
    %hippocampus
    %post
    [rho,p] = corr(results_pairs(:,2),hippo_cached_value(:,2), 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,hippo_cached_value(:,2), 'type', 'Spearman')

    %pre
    [rho,p] = corr(results_pairs(:,2),hippo_cached_value(:,1), 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,hippo_cached_value(:,1), 'type', 'Spearman')

    %pre-post
    diff_hippo_cached_value = hippo_cached_value(:,2) - hippo_cached_value(:,1);

    [rho,p] = corr(results_pairs(:,2),diff_hippo_cached_value, 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,diff_hippo_cached_value, 'type', 'Spearman')


    %vta
    %post
    [rho,p] = corr(results_pairs(:,2),vta_cached_value(:,2), 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,vta_cached_value(:,2), 'type', 'Spearman')

    %pre
    [rho,p] = corr(results_pairs(:,2),vta_cached_value(:,1), 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,vta_cached_value(:,1), 'type', 'Spearman')

    %pre-post
    diff_vta_cached_value = vta_cached_value(:,2) - vta_cached_value(:,1);

    [rho,p] = corr(results_pairs(:,2),diff_vta_cached_value, 'type', 'Spearman')
    [rho,p] = corr(diff_2A_2B,diff_vta_cached_value, 'type', 'Spearman')
    
    
    
    
end

