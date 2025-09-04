function T2 = T2_init(T2, FOV, system)

    %% T2 parameter initialization
    
    if ~isfield(T2, 'ref_phase')
        T2.ref_phase = 0;
    end

    if ~isfield(T2, 'exc_mode')
        T2.exc_mode = 'adiabatic_BIR4';
    end

    if strcmp(T2.exc_mode, 'adiabatic_AHP')
        if ~isfield(T2, 'ahp_exc_time')
            T2.ahp_exc_time = 3.0 *1e-3;
        end
        if ~isfield(T2, 'ahp_wmax')
            T2.ahp_wmax = 600 * 2*pi;
        end
    end    

    if strcmp(T2.exc_mode, 'adiabatic_BIR4')
        if ~isfield(T2, 'bir4_beta')
            T2.bir4_beta = 10;            
        end
        if ~isfield(T2, 'bir4_kappa')
            T2.bir4_kappa = atan(10);            
        end
        if ~isfield(T2, 'bir4_dw0')
            T2.bir4_dw0 = 30000;            
        end
        if ~isfield(T2, 'bir4_f1')
            T2.bir4_f1 = 640;            
        end
        if ~isfield(T2, 'bir4_alpha')
            T2.bir4_alpha = pi/2;            
        end
        if ~isfield(T2, 'bir4_dphi')
            T2.bir4_dphi = 0.175;            
        end
        if ~isfield(T2, 'bir4_tau')
            T2.bir4_tau = 0.01;            
        end
    end

    if ~isfield(T2, 'crush_nTwists')
        T2.crush_nTwists = 8;
    end
    
    T2.n_prep = numel(T2.prep_times);
    
    %% get T2 objects
    T2.T2_objs.d1                              = mr.makeDelay(system.gradRasterTime);    
    [T2.T2_objs.rf_90_td, T2.T2_objs.rf_90_tu] = T2_get_objs_EXC(T2, system);    
    T2.T2_objs.rf_comp_pos                     = T2_get_objs_COMP(system, T2.rfc_dur, T2.ref_phase + pi/2);
    T2.T2_objs.rf_comp_neg                     = T2_get_objs_COMP(system, T2.rfc_dur, T2.ref_phase - pi/2);    
    [T2.T2_objs.gx_crush, T2.T2_objs.gy_crush, T2.T2_objs.gz_crush] = T2_get_objs_CRUSH(T2, FOV, system);
    
    for j = 1:T2.n_prep
        temp_t_inter  = round( (T2.prep_times(j)-2*T2.rfc_dur) / 4 / system.gradRasterTime) * system.gradRasterTime;
        T2.T2_objs.t_inter(j) = mr.makeDelay(temp_t_inter);
        clear temp_t_inter;
    end
    
    %% calc total prep durations
    for j = 1 : T2.n_prep
        T2.duration(j) = T2_get_duration(T2, j);
    end

end

%% ----------------------------------------------------------------------------------------
function [rf_90_td, rf_90_tu] = T2_get_objs_EXC(T2, system)

    % AHP: tipdown and tipup
    if strcmp(T2.exc_mode, 'adiabatic_AHP')
        AHP      = AHP_get_params(T2.ahp_exc_time, T2.ahp_wmax);
        rf_90_td = AHP_pulse( system, 'tipdown', AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  T2.ref_phase + pi/2 );
        rf_90_tu = AHP_pulse( system, 'tipup',   AHP.tau, AHP.wmax, AHP.ratio, AHP.beta1, AHP.beta2,  T2.ref_phase + pi/2 );
    end
    
    % BIR4: tipdown and tipup
    if strcmp(T2.exc_mode, 'adiabatic_BIR4')
        rf_90_td = SIGPY_BIR4(T2.bir4_beta, T2.bir4_kappa, T2.bir4_dw0, T2.bir4_f1, T2.bir4_alpha, T2.bir4_dphi, T2.bir4_tau, T2.ref_phase + pi/2, 'tipdown', system);
        rf_90_tu = SIGPY_BIR4(T2.bir4_beta, T2.bir4_kappa, T2.bir4_dw0, T2.bir4_f1, T2.bir4_alpha, T2.bir4_dphi, T2.bir4_tau, T2.ref_phase + pi/2, 'tipup',   system);
    end

end

%% ----------------------------------------------------------------------------------------
function rf = T2_get_objs_COMP(system, tau, phaseOffset)

    % calc amplitude for refocusing pulse
    f1 = 1 / tau;
    
    % create block waveform for composite refocusing pulse
    N          = round(tau/system.rfRasterTime/4)*4;
    comp_amp   = ones(N,1) * f1;
    comp_phase = ones(N,1) * (-pi/2);
    comp_phase(N/4+1 : 3*N/4) = 0; % (-y +x -y) -> (+x)
    signal     = comp_amp .* exp(1i * comp_phase);
    
    % create pulseq rf struct from dummy
    rf             = mr.makeSincPulse( pi, system, 'Duration', N*system.rfRasterTime, 'timeBwProduct', 2, 'use', 'refocusing' );
    rf.signal      = signal;
    rf.phaseOffset = phaseOffset;

end

%% ----------------------------------------------------------------------------------------
function [gx_crush, gy_crush, gz_crush] = T2_get_objs_CRUSH(T2, FOV, system)

    % Crusher Gradients after preparation
    lim_grad       = 1/sqrt(3);  % limit for maximum gradient amplitude
    lim_slew       = 1/sqrt(3);  % limit for maximum slew rate
    
    % crush z
    crush_area_z   = T2.crush_nTwists / FOV.dz;  % [1/m]
    gz_crush       = mr.makeTrapezoid('z', 'Area', crush_area_z, 'maxGrad', system.maxGrad * lim_grad, 'maxSlew', system.maxSlew * lim_slew, 'system', system);
    
    % crush x
    crush_area_x   = T2.crush_nTwists / FOV.dx;  % [1/m]
    gx_crush       = mr.makeTrapezoid('x', 'Area', crush_area_x, 'maxGrad', system.maxGrad * lim_grad, 'maxSlew', system.maxSlew * lim_slew, 'system', system);
    gx_crush.delay = ceil(gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;
    
    % crush y
    crush_area_y   = T2.crush_nTwists / FOV.dy;  % [1/m]
    gy_crush       = mr.makeTrapezoid('y', 'Area', crush_area_y, 'maxGrad', system.maxGrad * lim_grad, 'maxSlew', system.maxSlew * lim_slew, 'system', system);
    gy_crush.delay = ceil(gz_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime + ...
                     ceil(gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;    

end

%% ----------------------------------------------------------------------------------------
function [duration] = T2_get_duration(T2, loop_T2)
duration =  mr.calcDuration(T2.T2_objs.rf_90_td) + ...
            mr.calcDuration(T2.T2_objs.t_inter(loop_T2)) + ...
            mr.calcDuration(T2.T2_objs.rf_comp_pos) + ...
            mr.calcDuration(T2.T2_objs.t_inter(loop_T2)) + ...
            mr.calcDuration(T2.T2_objs.t_inter(loop_T2)) + ...
            mr.calcDuration(T2.T2_objs.rf_comp_neg) + ...
            mr.calcDuration(T2.T2_objs.t_inter(loop_T2)) + ...
            mr.calcDuration(T2.T2_objs.rf_90_tu);
end