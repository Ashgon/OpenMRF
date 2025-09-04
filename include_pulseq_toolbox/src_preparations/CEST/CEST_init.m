function CEST = CEST_init(CEST, FOV, system)

%% calculate CEST pulse object and S0 reference delay
CEST.rf        = SL_get_SL_pulse(system, CEST.tau, CEST.f1, CEST.phase, 0);
CEST.ref_delay = mr.calcDuration(CEST.rf);

%% crusher gradients

lim_grad = 0.85;
lim_slew = 0.75;

CEST.crush_nTwist_x = 4;   % [] number of 2pi twists in read direction
CEST.crush_area_x   = CEST.crush_nTwist_x / FOV.dx;  % [1/m]
CEST.gx_crush       = mr.makeTrapezoid('x', 'Area', CEST.crush_area_x, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
CEST.tcrush_x       = mr.calcDuration(CEST.gx_crush);

CEST.crush_nTwist_y = 4;   % [] number of 2pi twists in read direction
CEST.crush_area_y   = CEST.crush_nTwist_y / FOV.dy;  % [1/m]
CEST.gy_crush       = mr.makeTrapezoid('y', 'Area', CEST.crush_area_y, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
CEST.tcrush_y       = mr.calcDuration(CEST.gy_crush);
CEST.gy_crush.delay = ceil(CEST.gx_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

CEST.crush_nTwist_z = 8;   % [] number of 2pi twists in slice thickness
CEST.crush_area_z   = CEST.crush_nTwist_z / FOV.dz;  % [1/m]
CEST.gz_crush       = mr.makeTrapezoid('z', 'Area', CEST.crush_area_z, 'maxGrad', system.maxGrad*lim_grad, 'maxSlew', system.maxSlew*lim_slew, 'system', system);
CEST.tcrush_z       = mr.calcDuration(CEST.gz_crush);    
CEST.gz_crush.delay = ceil(CEST.gy_crush.riseTime*1.05 / system.gradRasterTime) * system.gradRasterTime;

CEST.t_crush = max([mr.calcDuration(CEST.gx_crush); mr.calcDuration(CEST.gy_crush); mr.calcDuration(CEST.gz_crush)]);

end

