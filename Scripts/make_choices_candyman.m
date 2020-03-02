function [sim_choices] = make_choices_candyman(fit_params,version, mechanism, u_sub_r, r, rev_s, r_rev, choices, choice_trials, forcedchoices)
    
    sim_choices = NaN(1,120);
    
    if numel(fit_params) == 5
       alpha1 = fit_params(1);
       tau = fit_params(2);   
    elseif numel(fit_params)  == 6
       alpha1 = fit_params(1);
       alpha2 = fit_params(2);
       tau = fit_params(3);    
    elseif numel(fit_params)  == 7
       alpha1 = fit_params(1);
       alpha2 = fit_params(2);
       alpha3 = fit_params(3);
       tau = fit_params(4);
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
                    
                end
                
                
                %% revaluation - gives an estimate of values for devalued and revalued choice options (alpha2)
                inival_reval = values_phase1(:,rev_s);%both value vectors initialized with last value estimate from learning phase
                v_reval(1,:) = inival_reval;
                

                if version == 1 || version == 2 || version == 4
                    if numel(fit_params)  == 5
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha1*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(fit_params)  == 6
                           for t=1:numel(choices)
                               v_reval(t+1,:) = v_reval(t,:) + alpha2*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(fit_params)  == 7
                           %value updating in revaluation phase, differentiated for
                           %unchosen (alpha2) and chosen (alpha3) option
                           for t=1:numel(choices)
                               v_reval(t+1,r_rev(t,:)==-1) = v_reval(t,r_rev(t,:)==-1) + alpha2*(r_rev(t,r_rev(t,:)==-1) - v_reval(t,r_rev(t,:)==-1));
                               v_reval(t+1,r_rev(t,:)==1) = v_reval(t,r_rev(t,:)==1) + alpha3*(r_rev(t,r_rev(t,:)==1) - v_reval(t,r_rev(t,:)==1));
                           end
                    end
                
                else
                    if numel(fit_params)  == 5
                           for t=1:numel(choices)
                               nan_pos = isnan(r_rev(t,:));
                               v_reval(t+1,nan_pos) = v_reval(t,nan_pos);
                               v_reval(t+1,~nan_pos) = v_reval(t,~nan_pos) + alpha1*(r_rev(t,~nan_pos) - v_reval(t,~nan_pos));
                               %v_reval(t+1,:) = v_reval(t,:) + alpha1*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(fit_params)  == 6
                           for t=1:numel(choices)
                               nan_pos = isnan(r_rev(t,:));
                               v_reval(t+1,nan_pos) = v_reval(t,nan_pos);
                               v_reval(t+1,~nan_pos) = v_reval(t,~nan_pos) + alpha2*(r_rev(t,~nan_pos) - v_reval(t,~nan_pos));
                               %v_reval(t+1,:) = v_reval(t,:) + alpha2*(r_rev(t,:) - v_reval(t,:));
                           end
                    elseif numel(fit_params)  == 7
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
                

                %% forced choice phase
                %%use optimized parameters to generate values/value
                %%differences
                value_vector = values_phase1;
                value_vector(rev_s) = values_phase2;
  
                value_vector = value_vector(choice_trials);
                
                % turn values into softmax choice probs
                vd = value_vector(:,1)-value_vector(:,2); %value difference between left/right option
                p = 1./(1 + exp(-vd/tau)); %softmax rule
                pch = [p 1-p]; %choice probabilities for left/right option
                
                %compare softmax generated choice probabilities
                %(pch) with randomly drawn number and assign left or right
                %choice
                rand_choices=rand(length(pch),1);
                
                sim_choices(rand_choices<pch(:,1))=1;
                sim_choices(rand_choices>pch(:,1))=2;
                
                

end


