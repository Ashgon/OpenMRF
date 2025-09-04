function mri_signal = mg_sim_signals_fessler( Mxy, mask2D, cmaps, ktraj, NCoils, NR, NRead, Nxy )

% init
mri_signal = zeros(NR, NRead, NCoils);
resort_map = reshape(1:Nxy^2, Nxy,Nxy);
resort_map = resort_map(mask2D);
cmaps      = permute(cmaps, [2,3,1]);
ktraj      = ktraj / max(abs(ktraj(:)));
N          = Nxy * [1,1];

% simulate signals via fessler toolobx
parfor j=1:NR
    % images in real space
    temp_mxy   = Mxy(:,j);
    temp_image = zeros(Nxy,Nxy);
    temp_image(resort_map) = temp_mxy;
    temp_image = repmat(temp_image, 1, 1, NCoils) .* cmaps;

    % mri data in kspace -> fessler
    temp_ktraj  = squeeze(ktraj(:,j,:))';
    temp_signal = nufft(temp_image, nufft_init(temp_ktraj*pi, N, [6 6], N*2, N/2, 'minmax:kb'));
    % temp_signal = Gmri(temp_ktraj, true(N), 'fov', fov, 'basis', {'dirac'}, 'nufft', {N, [6 6], 2*N, N/2, 'table', 2^12, 'minmax:kb'}) * reshape(temp_image, [], NCoils); 

    % save signal of NR
    mri_signal(j,:,:) = temp_signal;
end

mri_signal = permute(mri_signal, [3,1,2]);

end

