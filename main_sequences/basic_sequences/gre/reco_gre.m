%%
clear

cmaps = [];
[Images, PULSEQ, study_info, cmaps, f0] = GRE_reco([], cmaps, [], []);

xtv(Images);
