function mri_signal = mg_sim_signals_phys(Mxy, T2s, df0, cmaps, mask1D, ktraj, fov, dwell, NCoils, NR, NRead, Nxy)

% define real space of phantom
[temp_x, temp_y] = meshgrid( -Nxy/2 : Nxy/2-1, -Nxy/2 : Nxy/2-1 );
temp_x = temp_x(:) * fov/Nxy;
temp_y = temp_y(:) * fov/Nxy;
temp_x(mask1D==0) = [];
temp_y(mask1D==0) = [];

% temp variables
temp_kx        = 2*pi * squeeze(ktraj(1,:,:));
temp_ky        = 2*pi * squeeze(ktraj(2,:,:));
temp_phi_off   = 2*pi * dwell * repmat(df0, 1, NR);
temp_t2s_decay = dwell ./ repmat(T2s, 1, NR);
temp_cmaps     = permute(repmat(cmaps, 1, 1, NR), [1, 3, 2]);

% case: NR=1
if NR==1
    temp_kx = temp_kx.';
    temp_ky = temp_ky.';
end

% init array for simulated mri signal
mri_signal = zeros(NRead, NR, NCoils);

% start simulation
tic
for j=1:NRead
        temp_signal = Mxy .* ...                            % simulated magnetization
                      exp(-temp_t2s_decay*j) .* ...         % T2s decay during adc
                      exp(-1i * ( ...
                      temp_phi_off*j + ...                  % dephasing: offresonance 
                      temp_x * temp_kx(:,j)' + ...          % dephasing: x-gradient
                      temp_y * temp_ky(:,j)' ));            % dephasing: y-gradient
        temp_signal = squeeze(sum( temp_cmaps .* repmat(temp_signal, 1, 1, NCoils) ));
        mri_signal(j,:,:) = temp_signal;
        disp([num2str(j/NRead*100,'%.2f') '%'])
end
toc

mri_signal = permute(mri_signal, [3, 2, 1]);

end

