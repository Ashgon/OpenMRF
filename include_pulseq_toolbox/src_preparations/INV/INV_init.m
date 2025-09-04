function INV = INV_init(INV, FOV, system)

%% magnetization inversion: global inversion pulse

    if ~isfield(INV,'rf_type')
        INV.rf_type = 'HYPSEC_inversion'; % default
    end    
   
    if strcmp(INV.rf_type, 'HYPSEC_inversion')
        if ~isfield(INV, 'tExc')
            INV.tExc = 10 *1e-3;  % [s] pulse duration
        end
        if ~isfield(INV, 'mu')
            INV.mu = 4.9;         % [ ] determines amplitude of frequency sweep
        end
        if ~isfield(INV, 'beta')
            INV.beta = 600;       % [Hz] peak amplitude
        end        
        [INV.rf] = SIGPY_HYPSEC(INV.beta, INV.mu, INV.tExc, 0, system);

    elseif strcmp(INV.rf_type, 'sinc')
        INV.tExc = 6 *1e-3;  % [s] pulse duration
        [INV.rf] = mr.makeSincPulse( pi, system, 'Duration', INV.tExc, 'SliceThickness', 1e6, 'use', 'inversion');

    elseif strcmp(INV.rf_type, 'sigpy_SLR')
        INV.tExc = 6 *1e-3;  % [s] pulse duration
        [INV.rf] = SIGPY_SLR(pi, INV.tExc, 0, 4, 'inv', 'ls', 0.01, 0.01, 0 ,1e6, system);

    elseif strcmp(INV.rf_type, 'AFP_inversion')
        ahp.tau      = 3 *1e-3;     % [s] half pulse duration
        ahp.wmax     = 600 * 2*pi;  % [rad/s]
        ahp          = AHP_get_params( ahp.tau, ahp.wmax);
        ahp.rf       = AFP_pulse( system, ahp.tau, ahp.wmax, ahp.ratio, ahp.beta1, ahp.beta2, 0 );
        INV.rf       = ahp.rf;
        INV.tExc     = 2*ahp.tau;
        clear ahp;
    end

    % crusher z
    INV.crush_nTwist = 8; % [] number of 2pi twists in read/phase/slice direction of voxel
    INV.crush_area_z = INV.crush_nTwist / FOV.dz;  % [1/m]
    INV.gz_crush     = mr.makeTrapezoid('z', 'Area', INV.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);

    % crusher x
    INV.crush_area_x   = INV.crush_nTwist / FOV.dx;  % [1/m]
    INV.gx_crush       = mr.makeTrapezoid('x', 'Area', INV.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);
    INV.gx_crush.delay = ceil(INV.gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

    % crusher y
    INV.crush_area_y   = INV.crush_nTwist / FOV.dy;  % [1/m]
    INV.gy_crush       = mr.makeTrapezoid('x', 'Area', INV.crush_area_z, 'maxGrad', system.maxGrad/sqrt(3), 'maxSlew', system.maxSlew/sqrt(3), 'system', system);
    INV.gy_crush.delay = ceil(INV.gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime + ...
                         ceil(INV.gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;
    INV.tcrush = mr.calcDuration(INV.gx_crush, INV.gy_crush, INV.gz_crush);

    % recovery delay
    if isfield(INV, 'inv_rec_time')
        if size(INV.inv_rec_time,2) > size(INV.inv_rec_time,1)
            INV.inv_rec_time = INV.inv_rec_time';
        end
        for j=1:numel(INV.inv_rec_time)
            INV.inv_rec_delay(j,1) = mr.makeDelay( round(INV.inv_rec_time(j)/system.gradRasterTime)*system.gradRasterTime );
        end
    else
        INV.inv_rec_time  = system.gradRasterTime;
        INV.inv_rec_delay = mr.makeDelay( round(INV.inv_rec_time/system.gradRasterTime)*system.gradRasterTime );
    end

end