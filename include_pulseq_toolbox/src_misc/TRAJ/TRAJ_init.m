function [TRAJ, FOV] = TRAJ_init(TRAJ, SPI, system)

%% set default values
if ~isfield(TRAJ, 'method')
    TRAJ.method = 'duyn';
end
if ~isfield(TRAJ, 'Nav')
    TRAJ.Nav = 10;
end
if ~isfield(TRAJ, 'Ndummy')
    TRAJ.Ndummy = 50;
end
if ~isfield(TRAJ, 'Trec')
    TRAJ.Trec = 150 *1e-3;
end
if ~isfield(TRAJ, 'slice_thickness')
    TRAJ.slice_thickness = 2 *1e-3;
end
if ~isfield(TRAJ, 'slice_offset')
    TRAJ.slice_offset = 50 *1e-3;
end
if ~isfield(TRAJ, 'exc_time')
    TRAJ.exc_time = 4 *1e-3;
end
if ~isfield(TRAJ, 'exc_tbw')
    TRAJ.exc_tbw = 4;
end
if ~isfield(TRAJ, 'exc_fa')
    TRAJ.exc_fa = 20 *pi/180;
end
if ~isfield(TRAJ, 'spoil_nTwists')
    TRAJ.spoil_nTwists = 8;
end
if ~isfield(TRAJ, 'lim_grad')
    TRAJ.lim_grad = 0.75;
end
if ~isfield(TRAJ, 'lim_slew')
    TRAJ.lim_slew = 0.75;
end
if ~isfield(TRAJ, 'mode_rewinder')
    TRAJ.mode_rewinder = 'off';
end

%% create adc for trajectory measurement
TRAJ.gx = SPI.gx;
TRAJ.gy = SPI.gy;
TRAJ.NR = numel(TRAJ.gx);

% only acquire during the spiral center-out gradients -> identical bandwidth & dwell time!
if strcmp(TRAJ.mode_rewinder, 'off')
    TRAJ.adc = SPI.adc;    
end

% also acquire the rewinder gradients -> this might change bandwidth & dwell time!
if strcmp(TRAJ.mode_rewinder, 'on')
    temp.geo.traj.t_adc     = mr.calcDuration(TRAJ.gx(1));
    temp.geo.traj.bandwidth = 1 / SPI.adc.dwell;
    temp.geo.traj.n_samples = temp.geo.traj.t_adc * temp.geo.traj.bandwidth;
    TRAJ.adc = SPI_calc_adc(temp, system);
    if TRAJ.adc.dwell ~= SPI.adc.dwell
        warning(['dwell time and bandwidth changed!']);
    end
    for j=1:TRAJ.NR % prevent timing errors
        TRAJ.gx(j).waveform  = [TRAJ.gx(j).waveform, zeros(1,100)];
        TRAJ.gy(j).waveform  = [TRAJ.gy(j).waveform, zeros(1,100)];
        TRAJ.gx(j).shape_dur = numel(TRAJ.gx(j).waveform)     * system.gradRasterTime;
        TRAJ.gy(j).shape_dur = numel(TRAJ.gy(j).waveform)     * system.gradRasterTime;
        TRAJ.gx(j).tt        = (1:numel(TRAJ.gx(j).waveform)) * system.gradRasterTime - system.gradRasterTime/2;
        TRAJ.gy(j).tt        = (1:numel(TRAJ.gy(j).waveform)) * system.gradRasterTime - system.gradRasterTime/2;
    end
end

TRAJ.adc.phaseOffset = pi/2;

%% calculate dummy delay object for reference scans
TRAJ.ref_delay = mr.makeDelay(mr.calcDuration(TRAJ.gx(1), TRAJ.gy(1)));

%% calculate objects for slice excitation
[TRAJ.rf, temp_exc] = mr.makeSincPulse( TRAJ.exc_fa, ...
                                                 system, ...
                                                 'Duration', TRAJ.exc_time,...
                                                 'SliceThickness', TRAJ.slice_thickness, ...
                                                 'maxSlew', system.maxSlew*TRAJ.lim_slew, ...
                                                 'timeBwProduct', TRAJ.exc_tbw, ...
                                                 'apodization', 0.5, ...
                                                 'PhaseOffset', 0, ...
										         'use', 'excitation' );
temp_reph = mr.makeTrapezoid('z', 'Area', -temp_exc.area/2, 'maxGrad', system.maxGrad*TRAJ.lim_grad, 'maxSlew', system.maxSlew*TRAJ.lim_slew, 'system', system);
TRAJ.rf.freqOffset =  TRAJ.slice_offset * temp_exc.amplitude;

%% set channels for gradients
for j=1:TRAJ.NR
    TRAJ.gx(j).channel = 'x';
    TRAJ.gy(j).channel = 'y';
end
TRAJ.g_exc_x  = temp_exc;  TRAJ.g_exc_x.channel  = 'x';
TRAJ.g_exc_y  = temp_exc;  TRAJ.g_exc_y.channel  = 'y';
TRAJ.g_reph_x = temp_reph; TRAJ.g_reph_x.channel = 'x';
TRAJ.g_reph_y = temp_reph; TRAJ.g_reph_y.channel = 'y';
clear temp_exc temp_reph;

%% calculate spoiler gradient
TRAJ.g_spoil_z = mr.makeTrapezoid('z', 'Area', TRAJ.spoil_nTwists / 1e-3, 'maxGrad', system.maxGrad*TRAJ.lim_grad, 'maxSlew', system.maxSlew*TRAJ.lim_slew, 'system', system);

%% recovery time between repetitions
TRAJ.Trec = mr.makeDelay(round(TRAJ.Trec/system.gradRasterTime)*system.gradRasterTime);

%% delay between objects
TRAJ.d1 = mr.makeDelay(system.gradRasterTime);

%% create dummy FOV object
FOV.Nxy      = 64;
FOV.fov_xy   = 2*TRAJ.slice_offset;
FOV.dz       = TRAJ.slice_thickness;
FOV.z_offset = 0;
FOV_init();

end

