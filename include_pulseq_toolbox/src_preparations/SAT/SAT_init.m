function SAT = SAT_init(SAT, FOV, system)

%% magnetization reset: global saturation pulse
if strcmp(SAT.mode, 'on')

    if ~isfield(SAT,'rf_type')
        SAT.rf_type = 'adiabatic_AHP'; % default
    end

    if strcmp(SAT.rf_type, 'adiabatic_AHP')
        if ~isfield(SAT, 'ahp')
            SAT.ahp.tau  = 3 *1e-3;     % [s] pulse duration
            SAT.ahp.wmax = 600 * 2*pi;  % [rad/s]
        end
        SAT.ahp  = AHP_get_params( SAT.ahp.tau, SAT.ahp.wmax);
        SAT.rf   = AHP_pulse( system, 'tipdown', SAT.ahp.tau, SAT.ahp.wmax, SAT.ahp.ratio, SAT.ahp.beta1, SAT.ahp.beta2, 0 );
        SAT.tExc = SAT.ahp.tau;

    elseif strcmp(SAT.rf_type, 'adiabatic_BIR4')
        if ~isfield(SAT, 'bir4_tau')
            SAT.bir4_tau = 10 *1e-3;            
        end
        if ~isfield(SAT, 'bir4_f1')
            SAT.bir4_f1 = 640;            
        end
        if ~isfield(SAT, 'bir4_beta')
            SAT.bir4_beta = 10;            
        end
        if ~isfield(SAT, 'bir4_kappa')
            SAT.bir4_kappa = atan(10);            
        end
        if ~isfield(SAT, 'bir4_dw0')
            SAT.bir4_dw0 = 30000;            
        end
        if ~isfield(SAT, 'bir4_alpha')
            SAT.bir4_alpha = pi/2;            
        end
        if ~isfield(SAT, 'bir4_dphi')
            SAT.bir4_dphi = 0.175;            
        end        
        SAT.rf = SIGPY_BIR4(SAT.bir4_beta, SAT.bir4_kappa, SAT.bir4_dw0, SAT.bir4_f1, SAT.bir4_alpha, SAT.bir4_dphi, SAT.bir4_tau, 0, 'tipdown', system);
    
    elseif strcmp(SAT.rf_type, 'sinc')
        SAT.tExc = 3 *1e-3;  % [s] pulse duration
        SAT.rf = mr.makeSincPulse( pi/2, system, 'Duration', SAT.tExc, 'use', 'saturation');

    elseif strcmp(SAT.rf_type, 'sigpy_SLR')
        SAT.tExc = 3 *1e-3;  % [s] pulse duration
        SAT.rf = SIGPY_SLR(pi/2, SAT.tExc, 0, 4, 'sat', 'ls', 0.01, 0.01, 0, 1e6, system); 

    end

    SAT.lim_grad = 1/sqrt(3);
    SAT.lim_slew = 1/sqrt(3);

    if ~isfield(SAT, 'crush_nTwist_x')
        SAT.crush_nTwist_x = 8;   % [] number of 2pi twists in x direction
    end
    if ~isfield(SAT, 'crush_nTwist_y')
        SAT.crush_nTwist_y = 8;   % [] number of 2pi twists in y direction
    end
    if ~isfield(SAT, 'crush_nTwist_z')
        SAT.crush_nTwist_z = 8;   % [] number of 2pi twists in z direction
    end

    SAT.crush_area_x   = SAT.crush_nTwist_x / FOV.dx;  % [1/m]
    SAT.gx_crush       = mr.makeTrapezoid('x', 'Area', SAT.crush_area_x, 'maxGrad', system.maxGrad*SAT.lim_grad, 'maxSlew', system.maxSlew*SAT.lim_slew, 'system', system);
    SAT.tcrush_x       = mr.calcDuration(SAT.gx_crush);

    SAT.crush_area_y   = SAT.crush_nTwist_y / FOV.dy;  % [1/m]
    SAT.gy_crush       = mr.makeTrapezoid('y', 'Area', SAT.crush_area_y, 'maxGrad', system.maxGrad*SAT.lim_grad, 'maxSlew', system.maxSlew*SAT.lim_slew, 'system', system);
    SAT.tcrush_y       = mr.calcDuration(SAT.gy_crush);

    SAT.crush_area_z   = SAT.crush_nTwist_z / FOV.dz;  % [1/m]
    SAT.gz_crush       = mr.makeTrapezoid('z', 'Area', SAT.crush_area_z, 'maxGrad', system.maxGrad*SAT.lim_grad, 'maxSlew', system.maxSlew*SAT.lim_slew, 'system', system);
    SAT.tcrush_z       = mr.calcDuration(SAT.gz_crush);    
    SAT.gz_crush.delay = ceil(SAT.gy_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;   

    % recovery delay
    if isfield(SAT, 'sat_rec_time')
        if size(SAT.sat_rec_time,2) > size(SAT.sat_rec_time,1)
            SAT.sat_rec_time = SAT.sat_rec_time';
        end
        for j=1:numel(SAT.sat_rec_time)
            SAT.sat_rec_delay(j,1) = mr.makeDelay( round(SAT.sat_rec_time(j)/system.gradRasterTime)*system.gradRasterTime );
        end
    else
        SAT.sat_rec_time  = system.gradRasterTime;
        SAT.sat_rec_delay = mr.makeDelay( round(SAT.sat_rec_time/system.gradRasterTime)*system.gradRasterTime );
    end

end

end