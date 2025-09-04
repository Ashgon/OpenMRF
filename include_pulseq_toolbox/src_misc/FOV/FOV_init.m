%% case: quadratic, isotropic
if ( ~isfield(FOV, 'Nx') && isfield(FOV, 'Nxy') )
    FOV.Nx = FOV.Nxy;
end
if ( ~isfield(FOV, 'Ny') && isfield(FOV, 'Nxy') )
    FOV.Ny = FOV.Nxy;
end
if ( ~isfield(FOV, 'fov_x') && isfield(FOV, 'fov_xy') )
    FOV.fov_x = FOV.fov_xy;
end
if ( ~isfield(FOV, 'fov_y') && isfield(FOV, 'fov_xy') )
    FOV.fov_y = FOV.fov_xy;
end

%% spatial resolution
if ~isfield(FOV, 'dx') 
    FOV.dx = FOV.fov_x / FOV.Nx;
end
if ~isfield(FOV, 'dy') 
    FOV.dy = FOV.fov_y / FOV.Nx;
end
if ~isfield(FOV, 'fov_z') 
    FOV.fov_z = FOV.dz;
end

%% standard orientation
% this might change if you rotate at the scanner!
FOV.orientation     = 'axial';
FOV.read_direction  = 'x';
FOV.phase_direction = 'y';
FOV.Rot_Mat         = [1 0 0; 0 1 0; 0 0 1];