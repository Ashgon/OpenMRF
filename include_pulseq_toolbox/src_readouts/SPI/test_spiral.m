%%
clear
pulseq_init();

%% FOV geometry
FOV.Nxy      = 128;         % [ ] matrix size
FOV.Nz       = 1;           % [ ] numer of "stack-of-spirals", 1 -> 2D
FOV.fov_xy   = 256  *1e-3;  % [m] FOV geometry
FOV.dz       = 5   *1e-3;   % [m] slab or slice thickness
FOV.z_offset = 0    *1e-3;  % [m] slice offset
FOV.fov_z    = FOV.dz;
FOV_init();

%% init spiral params
SPI_params_flash();
[SPI, ktraj_adc, ktraj_full, ktraj_reco] = SPI_init(SPI, FOV, system, 1);

%% add sequence blocks
for loop_kz = 1:FOV.Nz
for loop_NR = 1-SPI.Ndummy:SPI.NR
    SPI_add();
end
end

%% plot sequence and z trajecotry
seq.plot()

