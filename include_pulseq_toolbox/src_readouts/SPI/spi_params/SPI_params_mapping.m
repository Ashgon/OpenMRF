%% sequence params

% basic params
SPI.NR     = 8;             % [ ] number of repetitions
SPI.TR     = 0 *1e-3;       % [s] repetition time, 0 for minimization
SPI.Ndummy = 0;             % [ ] initial dummy loops

% slice excitation
SPI.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
SPI.exc_shape      = 'ex';          % only for sigpy: 'st' or 'ex' 
SPI.exc_time       = 2.0 *1e-3;     % [s] excitation time
SPI.exc_tbw        = 4;             % [ ] time bandwidth product
SPI.exc_fa_mode    = 'equal';       % 'equal',  'ramped',  'import'  
SPI.exc_flipangle  = 90 *pi/180;    % [rad] const FA or start of FA ramp -> set [] for auto mode
SPI.lim_gz_slew    = 0.8;           % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew  = 0.6;           % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 8;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 3.5 *1e-3;  % [s] time for spoiler and rewinder gradients
SPI.lim_spoil_slew = 0.8;        % [ ] reduce stimulation during gradient spoiling

% rf spoiling
SPI.spoil_rf_mode = 'quad';       % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 117 *pi/180;  % rf spoiling increment [rad]

%% spiral geoemtry mode
SPI.geo.traj_mode       = 'vds';         % 'standard' or 'vds'
SPI.geo.interleave_mode = 'Equal2Pi';    % 'Equal2Pi', 'Random2Pi', 'GoldenAngle', 'RandomGoldenAngle', 'RoundGoldenAngle'
SPI.deltak              = 1/FOV.fov_xy;  % [1/m] kspace sampling
SPI.kmax                = SPI.deltak * FOV.Nxy/2;

% spiral geometry: vds-hargreaves-toolbox
if strcmp(SPI.geo.traj_mode, 'vds')
    SPI.geo.Nvds     = 7.9;           % number of vds-spirals for sampling the kspce center
    SPI.geo.BW       = 400 *1e3;      % [Hz] bandwidth of spiral acquisition
    SPI.geo.Fcoeff   = [1  0];        % [1 0] for archimedean (equal density), [1 -0.99] for logarithmic (variable density)
    SPI.geo.grad_lim = 1/sqrt(3);     % limit of gradient field strength
    SPI.geo.slew_lim = 1/sqrt(3);     % limit of slew rate
    SPI.geo.kmax     = SPI.kmax;      % determines resolution
    SPI.geo.t_dwell  = 1/SPI.geo.BW;  % [s] dwell time for spiral acquisition
end
