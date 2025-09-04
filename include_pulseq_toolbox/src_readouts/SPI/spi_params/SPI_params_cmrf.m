%% sequence params

% import params from MRF struct
SPI.nr             = MRF.nr;    % [ ] number of readouts per heart beat
SPI.NR             = MRF.NR;    % [ ] number of repetitions
SPI.mrf_import.TRs = MRF.TRs;   % [s] repetition times
SPI.mrf_import.FAs = MRF.FAs;   % [rad] flip angles

% slice excitation
SPI.exc_mode       = 'sinc';        % 'sinc' or 'sigpy_SLR'
SPI.exc_shape      = 'ex';          % only for sigpy: 'st' or 'ex' 
SPI.exc_time       = 0.8 *1e-3;     % [s] excitation time
SPI.exc_tbw        = 2;             % [ ] time bandwidth product
SPI.exc_fa_mode    = 'import';    % 'equal',  'ramped',  'import' 
SPI.lim_gz_slew    = 0.9;         % [ ] reduce stimulation during slice excitation
SPI.lim_reph_slew  = 0.9;         % [ ] reduce stimulation during slice rephaser

% gradient spoiling
SPI.spoil_nTwist   = 4;          % [ ] number of 2pi twists in z-direction, 0 for balanced
SPI.spoil_duration = 0.8 *1e-3;  % [s] time for spoiler and rewinder gradients !!!! 1.0ms
SPI.lim_spoil_slew = 0.9;        % [ ] reduce stimulation during gradient spoiling

% rf spoiling
SPI.spoil_rf_mode = 'lin';      % rf spoiling mode: 'lin' or 'quad'
SPI.spoil_rf_inc  = 0 *pi/180;  % rf spoiling increment [rad]

%% spiral geoemtry mode
SPI.geo.traj_mode       = 'vds';              % 'standard' or 'vds'
SPI.geo.interleave_mode = 'RoundGoldenAngle'; % 'Equal2Pi', 'Random2Pi', 'GoldenAngle', 'RandomGoldenAngle', 'RoundGoldenAngle'
SPI.Nunique             = 48;                 % number of unique projection angles
SPI.deltak              = 1/FOV.fov_xy;  % [1/m] kspace sampling
SPI.kmax                = SPI.deltak * FOV.Nxy/2;

% spiral geometry: vds-hargreaves-toolbox
if strcmp(SPI.geo.traj_mode, 'vds')
    SPI.geo.Nvds     = 24;            % number of vds-spirals for sampling the kspce center
    SPI.geo.BW       = 500 *1e3;      % [Hz] bandwidth of spiral acquisition
    SPI.geo.Fcoeff   = [1  -0.5];     % [1 0] for archimedean (equal density), [1 -0.99] for logarithmic (variable density)
    SPI.geo.grad_lim = 1/sqrt(3);     % limit of gradient field strength
    SPI.geo.slew_lim = 1/sqrt(3);     % limit of slew rate
    SPI.geo.kmax     = SPI.kmax;      % determines resolution
    SPI.geo.t_dwell  = 1/SPI.geo.BW;  % [s] dwell time for spiral acquisition
end
