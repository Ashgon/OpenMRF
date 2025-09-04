% reconstruction of Siemens spin-echo (SE)
% use for T2 mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% read rawdata from .dat
study_path = 'Q:/data/Pulseq/Rawdata/tomgr/CimaX_Freiburg/250823_Open_MRF__NIST_T1_T2_T1p';
study_name = { 'meas_MID00383_FID31844_t2_se_cor_TE_10ms.dat',
               'meas_MID00384_FID31845_t2_se_cor_TE_15ms.dat',
               'meas_MID00385_FID31846_t2_se_cor_TE_20ms.dat',
               'meas_MID00386_FID31847_t2_se_cor_TE_30ms.dat',
               'meas_MID00387_FID31848_t2_se_cor_TE_40ms.dat',
               'meas_MID00388_FID31849_t2_se_cor_TE_60ms.dat',
               'meas_MID00389_FID31850_t2_se_cor_TE_85ms.dat',
               'meas_MID00390_FID31851_t2_se_cor_TE_120ms.dat',
               'meas_MID00391_FID31852_t2_se_cor_TE_175ms.dat',
               'meas_MID00392_FID31853_t2_se_cor_TE_250ms.dat'};

for j=1:numel(study_name)
    twix_obj = mapVBVD(fullfile(study_path, study_name{j}), 'ignoreSeg', 'removeOS');
    twix_obj = twix_obj{end};
    twix_obj = squeeze(twix_obj.image());
    rawdata(j,:,:,:) = permute(twix_obj,[2,3,1]);
end

% ifft
Images_coils = kspace2image(rawdata);

% cmaps: use the first image as the reference
[cmaps, ~, ~] = mg_espirit_cmaps(squeeze(Images_coils(1,:,:,:)), [], [], [], []);

% coil combined images
[NInv, NCoils, Ny, Nx] = size(Images_coils);
Images = zeros(NInv, Ny, Nx);
for j=1:NInv
    Images(j,:,:,:) = squeeze(sum(squeeze(Images_coils(j,:,:,:)) .* conj(cmaps)));
end

% zero interpolation filling
zero_params.onoff  = 1;
zero_params.radius = 0.5;
zero_params.factor = 2.0;
Images = mg_zero_filling(Images, zero_params);

% fit mask
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', NInv);

%% T2 mapping

% rotate images to real axis: use the first image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(1,:,:)),[NInv,1,1])));

TE = [10 15 20 30 40 60 85 120 175 250] *1e-3;
[ T2_Map, ~, R2_Map ] = mg_map_T12p( Images, TE, mask );

%% vis results
t2lims = [0 1000] *1e-3;
t2cmp  = get_cmp('T2', 1000, 1);

figure()
ax1 = subplot(1,2,1);
imagesc(T2_Map, t2lims); axis image; axis off; colormap(ax1, t2cmp); colorbar; title('T1 map [s]');
ax2 = subplot(1,2,2);
imagesc(R2_Map, [0.8 1]); axis image; axis off; colormap(ax2, turbo(1000)); colorbar; title('R2 map');
