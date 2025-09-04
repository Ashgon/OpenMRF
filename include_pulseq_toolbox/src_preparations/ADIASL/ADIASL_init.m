function ADIASL = ADIASL_init(ADIASL, FOV, system)

%% adiabatic spin-lock preparation: global inversion pulse
if strcmp(ADIASL.mode, 'on')

    % set default parameters
    if ~isfield(ADIASL, 'phi')
        ADIASL.phi = 0;
    end
    if ~isfield(ADIASL, 'tau')
        ADIASL.tau = 0.015;
    end
    if ~isfield(ADIASL, 'f1max')
        ADIASL.f1max = 500;
    end
    if ~isfield(ADIASL, 'beta')
        ADIASL.beta = 5;
    end
    if ~isfield(ADIASL, 'dwmax')
        ADIASL.dfmax = 200;
    end
    if ~isfield(ADIASL, 'phi_list')
        ADIASL.phi_list = repmat([0 1 1 0   1 0 0 1   0 0 1 1   1 1 0 0], 1, 100)';
    end

    % create dummy pulse objects
    ADIASL.rf = mr.makeSincPulse( pi/2, system, 'duration', ADIASL.tau, 'use', 'other');
    ADIASL.rf.signal = ADIASL.rf.signal * 0;

    % calculate waveform of hyperbolic sechans:
    ADIASL.f1       = ADIASL.f1max * sech(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1));
    ADIASL.df_down  = -2*ADIASL.dfmax*tanh(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1)); % !!!
    ADIASL.df_up    = -2*ADIASL.dfmax*tanh(ADIASL.beta*(2*ADIASL.rf.t/ADIASL.tau-1)); % improved performance if both pulses are identical; compare to Coletti et al!
    ADIASL.phi_down = cumsum(2*pi*ADIASL.df_down) * system.rfRasterTime;
    ADIASL.phi_up   = cumsum(2*pi*ADIASL.df_up)   * system.rfRasterTime;
    ADIASL.phi_down = wrapTo2Pi(ADIASL.phi_down - ADIASL.phi_down(round(numel(ADIASL.f1)/2)));
    ADIASL.phi_up   = wrapTo2Pi(ADIASL.phi_up   - ADIASL.phi_up(  round(numel(ADIASL.f1)/2)));

    % calc prep times
    ADIASL.adia_prep_times = ADIASL.N_HS * mr.calcDuration(ADIASL.rf);

    % crusher z
    ADIASL.crush_nTwist = 8; % [] number of 2pi twists in read/phase/slice direction of voxel
    ADIASL.crush_area_z = ADIASL.crush_nTwist / FOV.dz;  % [1/m]
    ADIASL.gz_crush     = mr.makeTrapezoid('z', 'Area', ADIASL.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);

    % crusher x
    ADIASL.crush_area_x   = ADIASL.crush_nTwist / FOV.dx;  % [1/m]
    ADIASL.gx_crush       = mr.makeTrapezoid('x', 'Area', ADIASL.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);
    ADIASL.gx_crush.delay = ceil(ADIASL.gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

    % crusher y
    ADIASL.crush_area_y   = ADIASL.crush_nTwist / FOV.dy;  % [1/m]
    ADIASL.gy_crush       = mr.makeTrapezoid('x', 'Area', ADIASL.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);
    ADIASL.gy_crush.delay = ceil(ADIASL.gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime + ...
                            ceil(ADIASL.gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

    ADIASL.tcrush = mr.calcDuration(ADIASL.gx_crush, ADIASL.gy_crush, ADIASL.gz_crush);
    
end

end