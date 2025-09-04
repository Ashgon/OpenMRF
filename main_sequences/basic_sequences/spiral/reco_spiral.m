%% reco 2d
clear
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco([], cmaps, [], [], []);
xtv(Images)

%% reco 3d
clear
cmaps = [];
[Images, PULSEQ, study_info, cmaps, dcf2D, f0, rawdata_noise] = SPI_reco_3D([], cmaps, [], [], []);
xtv(Images)