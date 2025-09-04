%% init pulseq
% basic SPI (spiral) readout
clear
seq_name = 'spi_single';

% optional flags
flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
flag_pns    = 0; % 0: off,  1: simulate PNS stimulation
flag_sound  = 0; % 0: off,  1: simulate gradient sound
flag_mrf    = 0; % 0: off,  1: simulate sequence via MRF toolbox

% init system, seq object and load pulseq user information
pulseq_init();

%% FOV geometry
FOV.Nxy      = 96;          % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% spiral sequence parameters
SPI_params_single_shot();
[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system);

SPI.Trec = 1;
SPI.Nrep = 5;

%% fat saturation and reset
FAT.mode = 'on';
FAT = FAT_init(FAT, FOV, system);

SAT.mode = 'on';
SAT = SAT_init(SAT, FOV, system);

%% noise pre-scans
SPI.Nnoise = 100;
SPI_add_prescans();

%% create sequence

% saturtion or crusher
SAT_add();   

% recovery time
seq.addBlock(mr.makeDelay(SPI.Trec));

for loop_rep = 1 : SPI.Nrep    

    % saturtion or crusher
    SAT_add();   
    
    % recovery time
    seq.addBlock(mr.makeDelay(SPI.Trec));

    % fat saturation
    FAT_add();

    % spiral imaging
    for loop_NR = 1-SPI.Ndummy : SPI.NR
        SPI_add();
    end
   
end

%% plot sequence diagram
seq.plot()

%% set definitions, check timings/gradients and export/backup files
filepath = [mfilename('fullpath') '.m'];
pulseq_exit();