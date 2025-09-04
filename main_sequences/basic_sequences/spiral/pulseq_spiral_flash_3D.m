%% init pulseq
% basic SPI (spiral) readout
clear
seq_name = 'spi_flash_3D';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 256;         % [ ] matrix size
FOV.Nz       = 48;          % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 96   *1e-3;  % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_flash();
SPI.NR             = 48;
SPI.spoil_nTwist   = 4 * FOV.Nz;
SPI.spoil_duration = 4.0 *1e-3;
SPI.kz_table       = 1:FOV.Nz;
[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

%% noise pre-scans
SPI.Nnoise = 100;
SPI_add_prescans();

%% create sequence
ndummy = SPI.Ndummy;

for loop_kz = SPI.kz_table
    for loop_NR = 1-ndummy : SPI.NR
        SPI_add();
    end
    ndummy = 0;
end

%% plot sequence diagram
% seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();