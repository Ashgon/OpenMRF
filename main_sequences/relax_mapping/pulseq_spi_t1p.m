%% init pulseq
% basis: SPI readout
% use for: T1rho, T2rho or T2 mapping
clear
seq_name = 't1p_mapping';

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
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_mapping();
[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

SPI.Trec = 5;

%% spin-lock params
SL.relax_type = {'T1p'};
SL.seq_type   = {'BSL'};
SL.tSL        = [8    16    24    32    40    48    56    64    72] *1e-3;
SL.fSL        = [350] * ones(numel(SL.tSL),1);
SL.inv_mode   = {'off'};

SL_params();
SL = SL_init(SL, FOV, system);

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% create sequence
for loop_SL = 1 : SL.nSL    
for loop_NR = 1-SPI.Ndummy : SPI.NR

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % sl preparation
    SL_add();

    % spiral imaging
    SPI_add();  

end    
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();