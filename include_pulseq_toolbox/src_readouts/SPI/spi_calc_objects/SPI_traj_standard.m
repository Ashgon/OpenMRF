function [traj] = SPI_traj_standard(kmax, t_adc, N, a, b, t_raster, n_samples)

    % function for calculating spiral trajecories
    % input parameters:
    % N: number of spiral loops, phi=0...N*2*pi
    % a: spiral geometry factor, 1: homogenous, <1: high pass, >1: low pass
    % b: spiral acquisition factor
    % t: time grid,  0.0 : system.gradRasterTime : adcTime;
    % output parameters: 2D kspace trajectory
    
    %% trajectory based on: Dipl. Kai Tobias Block
    % https://www.mpinat.mpg.de/597698/tblock_diploma.pdf
    n_raster = round(t_adc / t_raster);
    t        = linspace(0, t_adc, n_raster);
    phi      = 2*pi*N * sqrt(b+(1-b)*t_adc)/t_adc * t ./ sqrt(b+(1-b)*t);
    dphix    = pi/2 * exp(-phi*2);
    dphiy    = 0;
    kr       = kmax / (exp(2*pi*a*N)-1) * (exp(a*phi)-1); % log spiral
    kx       = kr .* cos( phi - dphix );
    ky       = kr .* sin( phi - dphiy );
    
    %% output vds trajectory
    traj.mode      = 'standard';
    traj.k         = kx + 1i * ky;
    traj.n_samples = n_samples;
    traj.t_adc     = t_adc;
    traj.t_dwell   = round( t_adc/n_samples / 100e-9 ) *100e-9;
    traj.bandwidth = 1/traj.t_dwell;

end

