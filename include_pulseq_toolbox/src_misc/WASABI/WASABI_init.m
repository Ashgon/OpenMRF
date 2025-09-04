function WASABI = WASABI_init(WASABI, FOV, system)

%% calculate WASABI pulse object and S0 reference delay
WASABI.rf = mr.makeBlockPulse( pi/10, system, ...
                               'Duration', WASABI.tau, ...
                               'freqOffset', 0, ...
                               'PhaseOffset', WASABI.phase, ...
                               'use', 'preparation' );
WASABI.rf.signal = WASABI.rf.signal*0 + WASABI.f1;
WASABI.ref_delay = mr.calcDuration(WASABI.rf);

%% crusher gradients

lim_grad = 1/sqrt(3);
lim_slew = 1/sqrt(3);

WASABI.crush_nTwist_x = 4;   % [] number of 2pi twists in read direction
WASABI.crush_area_x   = WASABI.crush_nTwist_x / FOV.dx;  % [1/m]
WASABI.gx_crush       = mr.makeTrapezoid('x', 'Area', WASABI.crush_area_x, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
WASABI.tcrush_x       = mr.calcDuration(WASABI.gx_crush);

WASABI.crush_nTwist_y = 4;   % [] number of 2pi twists in read direction
WASABI.crush_area_y   = WASABI.crush_nTwist_y / FOV.dy;  % [1/m]
WASABI.gy_crush       = mr.makeTrapezoid('y', 'Area', WASABI.crush_area_y, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
WASABI.tcrush_y       = mr.calcDuration(WASABI.gy_crush);
WASABI.gy_crush.delay = ceil(WASABI.gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

WASABI.crush_nTwist_z = 8;   % [] number of 2pi twists in slice thickness
WASABI.crush_area_z   = WASABI.crush_nTwist_z / FOV.dz;  % [1/m]
WASABI.gz_crush       = mr.makeTrapezoid('z', 'Area', WASABI.crush_area_z, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
WASABI.tcrush_z       = mr.calcDuration(WASABI.gz_crush);    
WASABI.gz_crush.delay = ceil(WASABI.gy_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

WASABI.t_crush = max([mr.calcDuration(WASABI.gx_crush); mr.calcDuration(WASABI.gy_crush); mr.calcDuration(WASABI.gz_crush)]);

end

