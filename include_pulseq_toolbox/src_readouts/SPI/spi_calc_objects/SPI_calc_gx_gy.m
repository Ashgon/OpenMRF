function [gx, gy, sx, sy] = SPI_calc_gx_gy(SPI, phi, system)

    % import k(t) from SPI.geo
    k    = SPI.geo.traj.k;
    k    = k * exp(-1i * atan2(imag(k(end)),real(k(end))));
    k(1) = 0;

    % combine k(t) with rewinder trajectory
    if ~strcmp(SPI.geo.traj_mode, 'import')
        k = k * SPI.kmax / abs(k(end));
        k = SPI_traj_loop_in(k, SPI.spoil_duration, 0.9, system);
    end
    
    % rotate k(t) around phi
    k = k * exp(1i * phi);

    % split k(t) to kx(t) and ky(t)
    kxy(1,:) = real(k);
    kxy(2,:) = imag(k);

    % calculate gradient waveforms and slew rates
    [gxy, sxy] = mr.traj2grad(kxy, 'RasterTime', system.gradRasterTime);
    
    % adc padding: start the gx gy gradients after adc dead time
    gxy = [zeros(2, 2 + round(SPI.adcTPad/system.gradRasterTime)), gxy];
    
    % save combined gradients: spiral + rewinder/spoiler
    gx       = mr.makeArbitraryGrad('x', gxy(1,:), 'first', gxy(1,1), 'last', gxy(1,end), 'system', system);
    gy       = mr.makeArbitraryGrad('y', gxy(2,:), 'first', gxy(2,1), 'last', gxy(2,end), 'system', system);
    gx.delay = system.gradRasterTime;
    gy.delay = system.gradRasterTime;
    gx.first = 0;
    gy.first = 0;
    gx.last  = 0;
    gy.last  = 0;

    % save slew rates
    sx = sxy(1,:)';
    sy = sxy(2,:)';
       
end

