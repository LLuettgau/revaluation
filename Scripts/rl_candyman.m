function [energy, chprob] = rl_candyman(version, mechanism, u_sub_r, r, rev_s, r_rev, r_reco, choices, choice_trials, forcedchoices, params)

if isnan(r_reco)
    if numel(params) == 2
       alpha1 = params(1);
       tau = params(2);   
    elseif numel(params)  == 3
       alpha1 = params(1);
       alpha2 = params(2);
       tau = params(3);    
    elseif numel(params)  == 4
       alpha1 = params(1);
       alpha2 = params(2);
       alpha3 = params(3);
       tau = params(4);
    end
else
     if numel(params) == 3
       alpha1 = params(1);
       alpha4 = params(2);
       tau = params(3);   
    elseif numel(params)  == 4
       alpha1 = params(1);
       alpha2 = params(2);
       alpha4 = params(3);
       tau = params(4);    
    elseif numel(params)  == 5
       alpha1 = params(1);
       alpha2 = params(2);
       alpha3 = params(3);
       alpha4 = params(4);
       tau = params(5);
    end
end    

   %% update values in learning phase - single updating with Rescorla-Wagner algorithm
      inival = 0.5; %value vectors initialized with .5
      v = zeros(size(r,1),6); %set up value matrix (rewards, stimuli)
      v(1,:) = inival;
                
                if mechanism == 1
                    % update values in learning phase
                    for t=1:size(r,1)
                        nan_pos = isnan(r(t,:));
                        v(t+1,nan_pos) = v(t,nan_pos);
                        v(t+1,~nan_pos) = v(t,~nan_pos) + alpha1*(r(t,~nan_pos) - v(t,~nan_pos));
                    end
                    %only take final values from learning phase to revaluation
                    values_phase1 = v(end,:);
                    
                elseif mechanism == 2
                    % update values in learning phase
                    for t=1:size(r,1)
                        v(t+1,:) = v(t,:) + alpha1*(r(t,:) - v(t,:));
                    end
                    %only take final values from learning phase to revaluation
                    values_phase1 = v(end,:);
                    values_phase1 = values_phase1.*u_sub_r;
                    
                elseif mechanism == 3
                    % update values in learning phase and scale per trial by reward
                    % subjective value
                    for t=1:size(r,1)
                        v(t+1,:) = v(t,:) + alpha1*(r(t,:) - v(t,:));
                        v(t+1,:) = v(t+1,:).*u_sub_r;
                    end
                    %only take final values from learning phase to revaluation
                    values_phase1 = v(end,:);  
                end
                
                
                %% revaluation - gives an estimate of values for devalued and revalued choice options (alpha2)
                inival_reval = values_phase1(:,rev_s);%both value vectors initialized with last value estimate from learning phase
                v_reval(1,:) = inival_reval;
                
                if isnan(r_reco)

                if version == 1 || version == 2 || version == 4
                    if numel(params)  == 2
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha1*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 3
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha2*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 4
                           %value updating in revaluation phase, differentiated for
                           %unchosen (alpha2) and chosen (alpha3) option
                           for t=1:numel(choices)
                               v_reval(t+1,r_rev(t,:)==-1) = v_reval(t,r_rev(t,:)==-1) + alpha2*(r_rev(t,r_rev(t,:)==-1) - v_reval(t,r_rev(t,:)==-1));
                               v_reval(t+1,r_rev(t,:)==1) = v_reval(t,r_rev(t,:)==1) + alpha3*(r_rev(t,r_rev(t,:)==1) - v_reval(t,r_rev(t,:)==1));
                           end
                    end
                
                else
                    if numel(params)  == 2
                           for t=1:numel(choices)
                               nan_pos = isnan(r_rev(t,:));
                               v_reval(t+1,nan_pos) = v_reval(t,nan_pos);
                               v_reval(t+1,~nan_pos) = v_reval(t,~nan_pos) + alpha1*(r_rev(t,~nan_pos) - v_reval(t,~nan_pos));
                               %v_reval(t+1,:) = v_reval(t,:) + alpha1*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 3
                           for t=1:numel(choices)
                               nan_pos = isnan(r_rev(t,:));
                               v_reval(t+1,nan_pos) = v_reval(t,nan_pos);
                               v_reval(t+1,~nan_pos) = v_reval(t,~nan_pos) + alpha2*(r_rev(t,~nan_pos) - v_reval(t,~nan_pos));
                               %v_reval(t+1,:) = v_reval(t,:) + alpha2*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 4
                           %value updating in revaluation phase, differentiated for
                           %unchosen (alpha2) and chosen (alpha3) option
                           for t=1:numel(choices)
                               %nan_pos(t,:) = isnan(r_rev(t,[1 3]));
                               %if t >= 2
                               %v_reval(t+1,nan_pos(t,:)) = v_reval(t,nan_pos(t,:));
                               %end
                               v_reval(t+1,isnan(r_rev(t,:))) = v_reval(t,isnan(r_rev(t,:)));      
                               v_reval(t+1,r_rev(t,:)==-1) = v_reval(t,r_rev(t,:)==-1) + alpha2*(r_rev(t,r_rev(t,:)==-1) - v_reval(t,r_rev(t,:)==-1));       
                               v_reval(t+1,r_rev(t,:)==1) = v_reval(t,r_rev(t,:)==1) + alpha3*(r_rev(t,r_rev(t,:)==1) - v_reval(t,r_rev(t,:)==1));
                           end
                    end
                end
                    
                    
                %only take final values from revaluation phase to forced
                %choice
                values_phase2 = v_reval(end,:);                
                
                %reconsolidation
                elseif ~isnan(r_reco) 
                    
                    if numel(params)  == 3
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha1*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 4
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha2*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(params)  == 5
                           %value updating in revaluation phase, differentiated for
                           %unchosen (alpha2) and chosen (alpha3) option
                           for t=1:numel(choices)
                               v_reval(t+1,r_rev(t,:)==-1) = v_reval(t,r_rev(t,:)==-1) + alpha2*(r_rev(t,r_rev(t,:)==-1) - v_reval(t,r_rev(t,:)==-1));
                               v_reval(t+1,r_rev(t,:)==1) = v_reval(t,r_rev(t,:)==1) + alpha3*(r_rev(t,r_rev(t,:)==1) - v_reval(t,r_rev(t,:)==1));
                           end
                    end
                    
                %only take final values from revaluation phase to forced
                %choice
                values_phase2 = v_reval(end,:);
                
                % updating of associative strength in RS post
                inival_reco = values_phase1;
                inival_reco(rev_s) = values_phase2;
                v_reco(1,:) = inival_reco;
                
                for rec=1:size(r_reco,1)
                    v_reco(rec+1,:) = v_reco(rec,:) + alpha4*(r_reco(rec,:) - v_reco(rec,:));
                end
                
                value_vector = v_reco(end,:);
                end

                %% forced choice phase
                %%minimize negative log likelihood of the parameters (alpha1, alpha2 and tau) 
                %%given the participant's choice data in forced choice phase
                if isnan(r_reco)
                   value_vector = values_phase1;
                   value_vector(rev_s) = values_phase2;
                end
                
                value_vector = value_vector(choice_trials);
                
                % turn values into softmax choice probs
                vd = value_vector(:,1)-value_vector(:,2); %value difference between left/right option
                p = 1./(1 + exp(-vd/tau)); %softmax rule
                pch = [p 1-p]; %choice probabilities for left/right option
                
                %compare subject's choices with softmax generated choice probabilities
                %(pch) in a new vector subpch --> how likely is the model's
                %(Softmax) prediction of a choice given the observed
                %choices of a participant
                
                subpch = zeros(size(forcedchoices,1),1);

                for k=1:size(subpch,1)
                    if forcedchoices(k) == 1 %if subject chose left, take left softmax generated choice probability
                       subpch(k) = pch(k,1);
                    elseif forcedchoices(k)==2 %if subject chose right, take right softmax generated choice probability
                       subpch(k) = pch(k,2);
                    end
                end
 
if isnan(r_reco) 
    if numel(params) == 2              
        if alpha1 < 0 || alpha1 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end                
    elseif numel(params) == 3               
        if alpha1 < 0 || alpha1 > 1 || alpha2 < 0 || alpha2 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end         
    elseif numel(params) == 4
        if alpha1 < 0 || alpha1 > 1 || alpha2 < 0 || alpha2 > 1 || alpha3 < 0 || alpha3 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end
    end
else
    if numel(params) == 3             
        if alpha1 < 0 || alpha1 > 1 || alpha4 < 0 || alpha4 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end                
    elseif numel(params) == 4               
        if alpha1 < 0 || alpha1 > 1 || alpha2 < 0 || alpha2 > 1 || alpha4 < 0 || alpha4 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end         
    elseif numel(params) == 5
        if alpha1 < 0 || alpha1 > 1 || alpha2 < 0 || alpha2 > 1 || alpha3 < 0 || alpha3 > 1 || alpha4 < 0 || alpha4 > 1 || tau <= 0 || tau > 100 
           energy = 9999999;
        else
           energy = -sum(log(subpch)); %calculate the negative log likelihood of choices
        end
        end
end
    
chprob = pch;


