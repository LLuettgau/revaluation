README for data and analysis code for Luettgau, Tempelmann, Kaiser & Jocham
Use revaluation/Scripts/analysis_script.m to rerun analyses of Experiment 1, 2, 3 and fmri.
Use revaluation/Scripts/analysis_script_exp5.m to rerun behavioral analyses of Experiment 5.
Use revaluation/Scripts/Luettgau_RL.m and rl_candyman.m to rerun computational modelling of Experiment 1, 2, 3 and fmri.
Use revaluation/Scripts/Simulation_RL.m and make_choices_candyman.m to rerun simulations (/sim_data/) based on computational parameters (/params/) of Experiment 1, 2, 3 and fmri.
Data is organized in MATLAB structure arrays, one for each subject.

Details for the behavioral data files format:

.infoP indicates subID

** .FOC (Pavlovian learning) columns, in candyman_fmri_pre in fmri**
1 Ð trial number
2 Ð target (0 = blue/left, 1 = red/right square around CS)
3 Ð response to target (1 = left, 2 = right)
4 Ð response time to target in sec
5 Ð CS texture number (as created by Psychtoolbox)
6 Ð CS pre rating score
7 Ð US texture number (as created by Psychtoolbox)
8 Ð US value (pre rating)
9 Ð CS onset (relative to start time)
10 Ð ISI (fixation cross) onset
11 Ð US onset
12 Ð ITI onset
13 Ð ITI length in sec
14 Ð US presented (1) or not (0)
15 Ð empty
16 Ð trial followed by a pause?
17-20 Ð computer times for events in 9-14
21 Ð pause onset (relative to start time)
22 Ð pause onset (computer time)
23 Ð pause end (relative to start time)
24 Ð pause duration 
25 Ð pause end (computer time)
26 Ð number of CS-US pair presented (1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B)

** .devaluation (revaluation choices) columns **
1 Ð trial number
2 Ð left CS number (1-6), 1 = lower valued CS, 2 = higher valued CS, 3-6 are lures
3 Ð right CS number (1-6), 1 = lower valued CS, 2 = higher valued CS
4 Ð left CS texture number (as created by Psychtoolbox)
5 Ð left CS pre rating score
6 Ð left CS associated US texture number
7 Ð left CS associated US value (pre rating)
8 Ð right CS texture number (as created by Psychtoolbox)
9 Ð right CS pre rating score
10 Ð right CS associated US texture number
11 Ð right CS associated US value (pre rating)
12 Ð response (1 = left chosen, 2 = right chosen) 
13 Ð response time (onset of choice to response) in sec
14 Ð 22 computer times in sec

**  .forcedchoicekanjis (decision probe) columns **
1 Ð trial number
2 Ð left CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 Ð right CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
4 Ð left CS texture number (as created by Psychtoolbox)
5 Ð left CS pre rating score
6 Ð left CS associated US texture number
7 Ð left CS associated US value (pre rating)
8 Ð right CS texture number (as created by Psychtoolbox)
9 Ð right CS pre rating score
10 Ð right CS associated US texture number
11 Ð right CS associated US value (pre rating)
12 Ð response (1 = left chosen, 2 = right chosen) 
13 Ð response time (onset of choice to response) in sec
14 Ð 22 computer times in sec


fMRI data files

**  .prerepsup (PRE fMRI run) columns **
1 Ð trial number
2 Ð CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 Ð US number (1-3) 1 = US-, 2 = US0, 3 = US+
4 Ð CS texture number (as created by Psychtoolbox)
5 Ð US texture number (as created by Psychtoolbox)
6 Ð CS onset (relative to start time)
7 Ð ISI (fixation cross) onset
8 Ð US onset
9 Ð ITI onset
10 Ð ITI length in sec
11 Ð attentional control task probe presented (1) or not (0)
12 Ð attentional control task probe onset
13 Ð order of probe (1 = right/wrong, 2 = wrong/right)
14 Ð correct (1)/incorrect (0) probe
15 Ð response (1 = left chosen, 2 = right chosen) 
16 Ð response time to probe in sec 
17 Ð response time relative to start time
18 Ð rewarded/correct response (1) or punished/incorrect response/time out (-1)
19 Ð onset of displaying of selected probe answer
20 Ð trial followed by a pause?
21 Ð pause onset
22 Ð pause end
23 Ð computer time of 6
24 Ð computer time of 7
25 Ð computer time of 8
26 Ð computer time of 9
27 Ð computer time of 12
28 Ð computer time of 15
29 Ð computer time of 16
30 Ð computer time of 19
31 Ð computer time of pause onset
32 Ð computer time of pause end

**  .postrepsup (POST fMRI run) columns **
1 Ð trial number
2 Ð CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 Ð US number (1-3) 1 = US-, 2 = US0, 3 = US+
4 Ð CS texture number (as created by Psychtoolbox)
5 Ð US texture number (as created by Psychtoolbox)
6 Ð CS onset (relative to start time)
7 Ð ISI (fixation cross) onset
8 Ð US onset
9 Ð ITI onset
10 Ð ITI length in sec
11 Ð attentional control task probe presented (1) or not (0)
12 Ð attentional control task probe onset
13 Ð order of probe
14 Ð correct (1)/incorrect (0) probe
15 Ð response (1 = left chosen, 2 = right chosen) 
16 Ð response time to probe in sec 
17 Ð response time relative to start time
18 Ð rewarded/correct response (1) or punished/incorrect response/time out (-1)
19 Ð onset of displaying of selected probe answer
20 Ð trial followed by a pause?
21 Ð pause onset
22 Ð pause end
23 Ð computer time of 6
24 Ð computer time of 7
25 Ð computer time of 8
26 Ð computer time of 9
27 Ð computer time of 12
28 Ð computer time of 15
29 Ð computer time of 16
30 Ð computer time of 19
31 Ð computer time of pause onset
32 Ð computer time of pause end

** .trigger (log file of fMRI triggers) **
1 Ð PRE run initial three volume onsets before start of experiment (computer time)
2 Ð PRE run initial synchronization triggers after every 15th trial (computer time)
3 Ð POST run initial three volume onsets before start of experiment (computer time)
4 Ð POST run initial synchronization triggers after every 15th trial (computer time)

fMRI data/extracted parameters

** /rsa_stats/ ** 
- contains .mat-files with extracted similarity structures (RSA) from right OFC (ending with _R_OFC_ROI.mat) and left hippocampus (ending with _L_Hippocampus_ROI.mat)
- PRE run starts with pre_
- POST run starts with post_
- PRE to POST difference starts with diff_
- _all_corr_ files contain all 66 pairwise correlation coefficients vectorized for n=42 subjects (for statistics)
- _sim_mat_ files contain group averaged squareform correlation coefficients (for illustration/plotting)

** /extracted_parameter_estimates/ **
- contains extracted parameters from several ROIs and several contrasts
