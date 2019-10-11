% simple reinforcement learner (Rescorla-Wagner update) with learning rate
% alpha and softmax temperature as free paramters
%Depending on options selected, will perform a grid search with 30 steps in each dimension to obtain an
%informed initial value for function optimization using fmincon, or 0/1 initial values or a particle
%swarm optimization
%Depending on selected mechanism, either a standard reinforcement learner
%or a associative strength model, multiplying the final values from the
%value estimation by the subjective value for each foodcue will be applied
%
%==============================================================================
function []=Luettgau_RL(version, mechanism, alphas)

% clear all; close all; clc

%version = 1; %indicate task version
pso = 0; %particle swarm optimization yes (1) or no (0)
%alphas = 2; %select number of learning rates, 1, 2, 3
inivals = 2; %set initial value for optimization: 1, 0, or 2 (based on grid search)
%mechanism = 1; %set value estimation algorithm, 1 = Rescorla-Wagner with 0, .5 and 1 as inputs, 2 = Rescorla-Wagner with 0/1 as inputs (depending on reward received or not), 
               %final value estimates are multiplied with individual rating
               %of unconditioned stimulus during initial foodcue rating
               %3 = trial-by-trial scaling
recons = 0;


data_path = '...' %set path to exp1, exp 2, exp3 or fmri here
res_dir = '...' %specify folder to save matrix with computational parameters

if version == 1
       rev_s = [3 5]; %exp1
       data_path = fullfile([data_path, 'exp1/']);
       datadir = data_path;
       vp = dir(fullfile([data_path, 'candymanLOG*']));
elseif version == 2
       rev_s = [1 3]; %exp2
       data_path = fullfile([data_path, 'exp2/']);
       datadir = data_path;
       vp = dir(fullfile([data_path, 'candyman2LOG*']));
elseif version == 3
       rev_s = [1 3 5]; %exp3
       data_path = fullfile([data_path, 'exp3/']);
       datadir = data_path;
       vp = dir(fullfile([data_path, 'candyman3LOG*']));
elseif version == 4
       rev_s = [3 5]; %fmri
       data_path = fullfile([data_path, 'fmri/']);
       datadir = data_path;
       vp = dir(fullfile([data_path, 'candyman_fmriLOG_9*']));
       recons = 1; %reconsolidation during RS post modelled?
end

vp = {vp.name};

nsteps = 15; %30 define steps within the parameter space/numer of parameter constellations to be tested for initial grid search (starting point)
space_alpha1 = exp(linspace(log(0.01), log(1), nsteps)); %set up learning rate alpha1 during the learning phase in log space
space_alpha2 = exp(linspace(log(0.00001), log(1), nsteps)); %set up learning rate alpha2 during revaluation (unchosen) in log space
space_alpha3 = exp(linspace(log(0.00001), log(1), nsteps)); %set up learning rate alpha3 during revaluation (chosen) in log space
space_alpha4 = exp(linspace(log(0.00001), log(1), nsteps)); %set up learning rate alpha4 for reconsolidation of S2A-O2 asso during post RS in log space
space_tau = exp(linspace(log(0.01), log(3), nsteps)); %set up softmax parameter tau in log space

for iteration=1:numel(vp)
    disp('Current iteration is:'); disp(iteration); drawnow
    
    
    %% preparations
    clear trials
    clear r
    clear r_new
    clear r_rev
    
    tic %timing for each participant
    data = importdata([datadir, vp{iteration}]);
    
    if version == 4
       predatadir = data_path;
       predata = dir(fullfile([data_path,'candyman_fmri_preLOG_9*']));
       predata = {predata.name};
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
    if recons == 1
       recodata = data.postrepsup;
    end
    
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
    
    if recons == 1
%        posreco = find(recodata(:,2) == 1 & recodata(:,3) == 1 | recodata(:,2) == 2 & recodata(:,3) == 1 | ... 
%                       recodata(:,2) == 3 & recodata(:,3) == 2 | recodata(:,2) == 4 & recodata(:,3) == 2 | ... 
%                       recodata(:,2) == 5 & recodata(:,3) == 3 | recodata(:,2) == 6 & recodata(:,3) == 3);
%        recodata = recodata(posreco,:);
%        
%        recodata = recodata(:,2);
       r_reco = repmat(u_sub_r,20,1); %received inputs (rewards)
    end
    
    if recons == 0
       r_reco = NaN;
    end
    
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
        
if recons == 0
    if pso == 0
        if alphas == 1
           if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
                for bb=1:numel(space_alpha2)
                    grid(aa,bb) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa),space_tau(bb)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                end
            end
            
            gridopt = min(min(grid)); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
    
            opt_alpha1 = space_alpha1(d1);
            opt_tau = space_tau(d2);
    
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_tau;
            end
    
    
        elseif alphas == 2
            if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
%                 disp(aa); drawnow
                for bb=1:numel(space_alpha2)
                    for cc=1:numel(space_tau)
                        grid(aa,bb,cc) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa), space_alpha2(bb), space_tau(cc)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                    end
                end
            end
            
            gridopt = min(min(min(grid))); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2, d3] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
    
            opt_alpha1 = space_alpha1(d1);
            opt_alpha2 = space_alpha2(d2);
            opt_tau = space_tau(d3);
    
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_alpha2;
            comp_p.params(iteration,3) = opt_tau;
            end
        
        elseif alphas == 3
            if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps, nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
%                 disp(aa); drawnow
                for bb=1:numel(space_alpha2)
                    for cc=1:numel(space_alpha3)  
                        for dd=1:numel(space_tau)
                            grid(aa,bb,cc,dd) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa), space_alpha2(bb), space_alpha3(cc), space_tau(dd)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                        end
                    end
                end
            end

            gridopt = min(min(min(min(grid)))); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2, d3, d4] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
    
            opt_alpha1 = space_alpha1(d1);
            opt_alpha2 = space_alpha2(d2);
            opt_alpha3 = space_alpha3(d3);
            opt_tau = space_tau(d4);
    
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_alpha2;
            comp_p.params(iteration,3) = opt_alpha3;
            comp_p.params(iteration,4) = opt_tau;
            end
        end
            %function optimization
            x = zeros(1,alphas+1);
            lb = zeros(1,numel(x));
            ub = [ones(1,numel(x)-1), 100];
            
            if alphas == 1            
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2)]);
            elseif alphas == 2
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3)]);
            elseif alphas == 3
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4)]);
            end
            %tryout different initial values, 0, 1, or based on 30x30x30 grid
            %search
            if inivals == 2
               if alphas == 1
                  ini = [comp_p.params(iteration,1) comp_p.params(iteration,2)];
               elseif alphas == 2
                      ini = [comp_p.params(iteration,1) comp_p.params(iteration,2) comp_p.params(iteration,3)];
               elseif alphas == 3 
                      ini = [comp_p.params(iteration,1) comp_p.params(iteration,2) comp_p.params(iteration,3) comp_p.params(iteration,4)];
               end
                   
            elseif inivals == 1
               ini = ones(1,numel(x));
            elseif inivals == 0
               ini = zeros(1,numel(x));
            end
    
            [minimum,energy,EXITFLAG] = fmincon(fun, ini,[],[],[],[],lb,ub);%, options)
        
    elseif pso == 1
    %function optimization with pso
    x = zeros(1,alphas+1);
    lb = zeros(1,numel(x));
    ub = [ones(1,numel(x)-1), 100];
    
    options = optimoptions('particleswarm','SwarmSize',10000000);
    if alphas == 1            
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2)]);
    elseif alphas == 2
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3)]);
    elseif alphas == 3
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4)]);
    end
    
        [minimum,energy,EXITFLAG] = particleswarm(fun, numel(x), lb,ub);%, options)
    end
    
    cp.params(iteration,1:numel(x)) = minimum;
    cp.params(iteration,numel(x)+1) = energy;
    cp.params(iteration,numel(x)+2) = exp(-energy/size(choice_trials,1));
    cp.params(iteration,numel(x)+3) = size(choice_trials,1);

%     cp.params(iteration,numel(x)+3) = nansum(r(:,rev_s(1)));
%     cp.params(iteration,numel(x)+4) = nansum(r(:,rev_s(1)+1));
%     cp.params(iteration,numel(x)+5) = nansum(r(:,rev_s(1)+1))-nansum(r(:,rev_s(1)));
%     cp.params(iteration,numel(x)+6) = nansum(r(:,rev_s(2)));
%     cp.params(iteration,numel(x)+7) = nansum(r(:,rev_s(2)+1));
%     cp.params(iteration,numel(x)+8) = nansum(r(:,rev_s(2)+1))-nansum(r(:,rev_s(2)));
%     if version == 3
%        cp.params(iteration,numel(x)+9) = nansum(r(:,rev_s(3)));
%        cp.params(iteration,numel(x)+10) = nansum(r(:,rev_s(3)+1));
%        cp.params(iteration,numel(x)+11) = nansum(r(:,rev_s(3)+1))-nansum(r(:,rev_s(3))); 
%     end

%version 4 models (reconsolidation)
elseif recons == 1
    
     if pso == 0
        if alphas == 1
           if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
                for bb=1:numel(space_alpha4)
                    for cc=1:numel(space_tau)
                        grid(aa,bb,cc) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa),space_alpha4(bb),space_tau(cc)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                    end
                end
            end
            
            gridopt = min(min(min(grid))); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2, d3] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
            
                
            opt_alpha1 = space_alpha1(d1);
            opt_alpha4 = space_alpha4(d2);
            opt_tau = space_tau(d3);
    
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_alpha4;
            comp_p.params(iteration,3) = opt_tau;
            end
    
    
        elseif alphas == 2
            if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
                for bb=1:numel(space_alpha2)
                    for cc=1:numel(space_alpha4)  
                        for dd=1:numel(space_tau)
                            grid(aa,bb,cc,dd) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa), space_alpha2(bb), space_alpha4(cc), space_tau(dd)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                        end
                    end
                end
            end
            
            gridopt = min(min(min(min(grid)))); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2, d3, d4] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
    
            opt_alpha1 = space_alpha1(d1);
            opt_alpha2 = space_alpha2(d2);
            opt_alpha4 = space_alpha4(d3);
            opt_tau = space_tau(d4);
            
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_alpha2;
            comp_p.params(iteration,3) = opt_alpha4;
            comp_p.params(iteration,4) = opt_tau;
            end
        
        elseif alphas == 3
            if inivals == 2 %only perform grid search if grid search informed in inital values are needed
            %grid search to obtain initial values
            grid = zeros(nsteps, nsteps, nsteps, nsteps); 
            for aa = 1:numel(space_alpha1)
                for bb=1:numel(space_alpha2)
                    for cc=1:numel(space_alpha3)  
                        for dd=1:numel(space_alpha4)
                            for ee=1:numel(space_tau)
                                grid(aa,bb,cc,dd,ee) =  rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, [space_alpha1(aa), space_alpha2(bb), space_alpha3(cc), space_alpha3(dd), space_tau(ee)]); %write returned energy parameter of the RL function to grid of alpha/tau  
                            end
                        end
                    end
                end
            end

            gridopt = min(min(min(min(min(grid))))); %find the grid option (minimum), minimum negative log likelihood of the parameter
   
            gridopt_pos = find(grid==gridopt);
    
            [d1, d2, d3, d4, d5] = ind2sub(size(grid),  gridopt_pos(1)); %find coordinates of grid optimum (alpha1, alpha2, alpha3 & tau)
    
            opt_alpha1 = space_alpha1(d1);
            opt_alpha2 = space_alpha2(d2);
            opt_alpha3 = space_alpha3(d3);
            opt_alpha4 = space_alpha3(d4);
            opt_tau = space_tau(d5);
    
            comp_p.params(iteration,1) = opt_alpha1;
            comp_p.params(iteration,2) = opt_alpha2;
            comp_p.params(iteration,3) = opt_alpha3;
            comp_p.params(iteration,4) = opt_alpha4;
            comp_p.params(iteration,5) = opt_tau;
            end
        end
            %function optimization
            x = zeros(1,alphas+2);
            lb = zeros(1,numel(x));
            ub = [ones(1,numel(x)-1), 100];
            
            if alphas == 1            
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3)]);
            elseif alphas == 2
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4)]);
            elseif alphas == 3
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4) x(5)]);
            end
            %tryout different initial values, 0, 1, or based on 30x30x30 grid
            %search
            if inivals == 2
               if alphas == 1
                  ini = [comp_p.params(iteration,1) comp_p.params(iteration,2) comp_p.params(iteration,3)];
               elseif alphas == 2
                      ini = [comp_p.params(iteration,1) comp_p.params(iteration,2) comp_p.params(iteration,3) comp_p.params(iteration,4)];
               elseif alphas == 3 
                      ini = [comp_p.params(iteration,1) comp_p.params(iteration,2) comp_p.params(iteration,3) comp_p.params(iteration,4) comp_p.params(iteration,5)];
               end
                   
            elseif inivals == 1
               ini = ones(1,numel(x));
            elseif inivals == 0
               ini = zeros(1,numel(x));
            end
    
            [minimum,energy,EXITFLAG] = fmincon(fun, ini,[],[],[],[],lb,ub);%, options)
        
    elseif pso == 1
    %function optimization with pso
    x = zeros(1,alphas+2);
    lb = zeros(1,numel(x));
    ub = [ones(1,numel(x)-1), 100];
    
    options = optimoptions('particleswarm','SwarmSize',10000000);
    if alphas == 1            
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3)]);
    elseif alphas == 2
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4)]);
    elseif alphas == 3
               fun = @(x) rl_candyman(version,mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices,[x(1) x(2) x(3) x(4) x(5)]);
    end
    
        [minimum,energy,EXITFLAG] = particleswarm(fun, numel(x), lb,ub);%, options)
    end
    
    cp.params(iteration,1:numel(x)) = minimum;
    cp.params(iteration,numel(x)+1) = energy;
    cp.params(iteration,numel(x)+2) = exp(-energy/size(choice_trials,1));
    cp.params(iteration,numel(x)+3) = size(choice_trials,1);
    
    
    toc
end
end
end
    %cp.params(:,numel(x)+3:end) = [];
    pos_cp = find(cp.params(:,end)>0);
    cp.params = cp.params(pos_cp,:);
    
if version == 1
    cd(res_dir)
	if mechanism == 1
    	if pso == 0
           save([pwd,filesep,'exp1_cp_rw_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp1_cp_rw_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 2
        if pso == 0
           save([pwd,filesep,'exp1_cp_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp1_cp_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 3
        if pso == 0
           save([pwd,filesep,'exp1_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp1_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end   
    end
elseif version == 2
    cd(res_dir)
    if mechanism == 1
    	if pso == 0
           save([pwd,filesep,'exp2_cp_rw_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp2_cp_rw_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 2
        if pso == 0
           save([pwd,filesep,'exp2_cp_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp2_cp_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 3    
        if pso == 0
           save([pwd,filesep,'exp2_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp2_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    end
elseif version == 3
    cd(res_dir) 
    if mechanism == 1
    	if pso == 0
           save([pwd,filesep,'exp3_cp_rw_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp3_cp_rw_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 2
        if pso == 0
           save([pwd,filesep,'exp3_cp_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp3_cp_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 3    
        if pso == 0
           save([pwd,filesep,'exp3_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'exp1_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end 
        
    end   
elseif version == 4 && recons == 0
    cd(res_dir)
    if mechanism == 1
    	if pso == 0
           save([pwd,filesep,'fmri_cp_rw_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_rw_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 2
        if pso == 0
           save([pwd,filesep,'fmri_cp_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 3    
        if pso == 0
           save([pwd,filesep,'fmri_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_tbt_asso_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end 
        
    end  
    
elseif version == 4 && recons == 1
    cd(res_dir)
    if mechanism == 1
    	if pso == 0
           save([pwd,filesep,'fmri_cp_rw_recons_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_rw_recons_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 2
        if pso == 0
           save([pwd,filesep,'fmri_cp_asso_recons_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_asso_recons_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end
    elseif mechanism == 3    
        if pso == 0
           save([pwd,filesep,'fmri_cp_tbt_asso_recons_', num2str(alphas), 'Alpha_1Tau_', num2str(inivals), '_ini.mat'],'cp');
        else
           save([pwd,filesep,'fmri_cp_tbt_asso_recons_', num2str(alphas), 'Alpha_1Tau_pso.mat'],'cp');
        end 
        
    end  
end

end
