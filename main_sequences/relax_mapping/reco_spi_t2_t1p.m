% reconstruction of pulseq spiral mapping sequences
% use for T2/T1p mapping of the NIST phantom
% mgram; V1; 28.11.2024

%% start reco
clear

% reconstruction of spiral datasets
[Images, PULSEQ] = SPI_reco();

% fit mask
Nimages = size(Images,1);
[mask, mask3D] = mg_get_mask_fit(squeeze(mean(abs(Images))), 'holes', Nimages);

%% T2 or T1p mapping

% rotate images to real axis: use the first image as the reference
Images = real(Images .* exp(-1i*repmat(angle(Images(1,:,:)),[Nimages,1,1])));

if isfield(PULSEQ, 'T2')
    TE = PULSEQ.T2.prep_times;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, TE, mask );
end
if isfield(PULSEQ, 'SL')
    tSL = PULSEQ.SL.tSL;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, tSL, mask );
end
if isfield(PULSEQ, 'ADIASL')
    tSL = PULSEQ.ADIASL.adia_prep_times;
    [ T12p_Map, ~, R2_Map ] = mg_map_T12p( Images, tSL, mask );
end

% vis results
figure()
subplot(1,2,1)
imagesc(T12p_Map, [0 1]); axis image; axis off; colormap(turbo(1000));
subplot(1,2,2)
imagesc(R2_Map, [0.8 1]); axis image; axis off; colormap(turbo(1000));