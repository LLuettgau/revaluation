# revaluation
Data and code underlying manuscript "Decisions bias future choices by modifying hippocampal associative memories" by Lennart Luettgau, Claus Tempelmann, Luca Franziska Kaiser, Gerhard Jocham

Use revaluation/Scripts/analysis_script.m to rerun analyses of Experiment 1, 2, 3 and fmri.

Use revaluation/Scripts/analysis_script_exp4.m to rerun behavioral analyses of Experiment 4.

Use revaluation/Scripts/Luettgau_RL.m and rl_candyman.m to rerun computational modelling of Experiment 1, 2, 3 and fmri.

Use revaluation/Scripts/Simulation_RL.m and make_choices_candyman.m to rerun simulations (/sim_data/) based on computational parameters (/params/) of Experiment 1, 2, 3 and fmri.

Data is organized in MATLAB structure arrays, one for each subject.

Details for the behavioral data files format:

.infoP indicates subID

** .FOC (Pavlovian learning) columns, in candyman_fmri_pre in fmri**
1 – trial number
2 – target (0 = blue/left, 1 = red/right square around CS)
3 – response to target (1 = left, 2 = right)
4 – response time to target in sec
5 – CS texture number (as created by Psychtoolbox)
6 – CS pre rating score
7 – US texture number (as created by Psychtoolbox)
8 – US value (pre rating)
9 – CS onset (relative to start time)
10 – ISI (fixation cross) onset
11 – US onset
12 – ITI onset
13 – ITI length in sec
14 – US presented (1) or not (0)
15 – empty
16 – trial followed by a pause?
17-20 – computer times for events in 9-14
21 – pause onset (relative to start time)
22 – pause onset (computer time)
23 – pause end (relative to start time)
24 – pause duration 
25 – pause end (computer time)
26 – number of CS-US pair presented (1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B)

** .devaluation (revaluation choices) columns **
1 – trial number
2 – left CS number (1-6), 1 = lower valued CS, 2 = higher valued CS, 3-6 are lures
3 – right CS number (1-6), 1 = lower valued CS, 2 = higher valued CS
4 – left CS texture number (as created by Psychtoolbox)
5 – left CS pre rating score
6 – left CS associated US texture number
7 – left CS associated US value (pre rating)
8 – right CS texture number (as created by Psychtoolbox)
9 – right CS pre rating score
10 – right CS associated US texture number
11 – right CS associated US value (pre rating)
12 – response (1 = left chosen, 2 = right chosen) 
13 – response time (onset of choice to response) in sec
14 – 22 computer times in sec

**  .forcedchoicekanjis (decision probe) columns **
1 – trial number
2 – left CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 – right CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
4 – left CS texture number (as created by Psychtoolbox)
5 – left CS pre rating score
6 – left CS associated US texture number
7 – left CS associated US value (pre rating)
8 – right CS texture number (as created by Psychtoolbox)
9 – right CS pre rating score
10 – right CS associated US texture number
11 – right CS associated US value (pre rating)
12 – response (1 = left chosen, 2 = right chosen) 
13 – response time (onset of choice to response) in sec
14 – 22 computer times in sec


fMRI data files

**  .prerepsup (PRE fMRI run) columns **
1 – trial number
2 – CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 – US number (1-3) 1 = US-, 2 = US0, 3 = US+
4 – CS texture number (as created by Psychtoolbox)
5 – US texture number (as created by Psychtoolbox)
6 – CS onset (relative to start time)
7 – ISI (fixation cross) onset
8 – US onset
9 – ITI onset
10 – ITI length in sec
11 – attentional control task probe presented (1) or not (0)
12 – attentional control task probe onset
13 – order of probe (1 = right/wrong, 2 = wrong/right)
14 – correct (1)/incorrect (0) probe
15 – response (1 = left chosen, 2 = right chosen) 
16 – response time to probe in sec 
17 – response time relative to start time
18 – rewarded/correct response (1) or punished/incorrect response/time out (-1)
19 – onset of displaying of selected probe answer
20 – trial followed by a pause?
21 – pause onset
22 – pause end
23 – computer time of 6
24 – computer time of 7
25 – computer time of 8
26 – computer time of 9
27 – computer time of 12
28 – computer time of 15
29 – computer time of 16
30 – computer time of 19
31 – computer time of pause onset
32 – computer time of pause end

**  .postrepsup (POST fMRI run) columns **
1 – trial number
2 – CS number (1-6), 1 = CS-A, 2= CS-B, 3 = CS0A, 4 = CS0B, 5 = CS+A, 6 = CS+B
3 – US number (1-3) 1 = US-, 2 = US0, 3 = US+
4 – CS texture number (as created by Psychtoolbox)
5 – US texture number (as created by Psychtoolbox)
6 – CS onset (relative to start time)
7 – ISI (fixation cross) onset
8 – US onset
9 – ITI onset
10 – ITI length in sec
11 – attentional control task probe presented (1) or not (0)
12 – attentional control task probe onset
13 – order of probe
14 – correct (1)/incorrect (0) probe
15 – response (1 = left chosen, 2 = right chosen) 
16 – response time to probe in sec 
17 – response time relative to start time
18 – rewarded/correct response (1) or punished/incorrect response/time out (-1)
19 – onset of displaying of selected probe answer
20 – trial followed by a pause?
21 – pause onset
22 – pause end
23 – computer time of 6
24 – computer time of 7
25 – computer time of 8
26 – computer time of 9
27 – computer time of 12
28 – computer time of 15
29 – computer time of 16
30 – computer time of 19
31 – computer time of pause onset
32 – computer time of pause end

** .trigger (log file of fMRI triggers) **
1 – PRE run initial three volume onsets before start of experiment (computer time)
2 – PRE run initial synchronization triggers after every 15th trial (computer time)
3 – POST run initial three volume onsets before start of experiment (computer time)
4 – POST run initial synchronization triggers after every 15th trial (computer time)

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
