%% install OpenMRF packages
disp(' ');
disp('     start installing packages for OpenMRF ...');
disp(' ');

%% include packages
addpath(genpath('include_pulseq_original'));
disp('     ... original pulseq repository included!');

addpath(genpath('include_pulseq_toolbox'));
disp('     ... pulseq toolbox functions included!');

addpath(genpath('user_specifications'));
disp('     ... pulseq user specifications included!');

addpath(genpath('include_pre_sim_library'));
disp('     ... pre-sim library included!');

addpath(genpath('include_misc'));
disp('     ... misc packages included!');

addpath(genpath('include_miitt/src_lr_reco'));
disp('     ... miitt low rank reco included!');

disp('     ... install michigan image reconstruction toolbox (MIRT)');
run('include_miitt/src_mirt/setup.m');
disp(' ');

%% check pulseq information
[pulseq_user, pulseq_lab, pulseq_path, pulseq_scanner] = pulseq_get_user_definitions([], 1);
disp(' ');
disp('     ... have fun!');
disp(' ');

%% check and display duplicate functions
pulseq_find_duplicates(mfilename('fullpath'));
clear;