# mgram 11.07.2025

# Contens:
functions for generating SLR or adiabatic pulses using the sigpy toolbox.
step 1: install sigpy to a local python environment (https://github.com/mikgroup/sigpy/tree/main)
step 2: check your new sigpy environment path and the python command, e.g.:
windows: env -> 'C:\\ProgramData\\miniforge3\\Scripts\\activate.bat SIGPY'
windows: python -> 'C:\\Users\\...\\.conda\\envs\\SIGPY\\python'
unix: env -> '/home/.../miniforge3/bin/activate SIGPY'
unix: python -> 'python3'
step 3: copy your env path and python command to the init_SIGPY_windows_unix.m script
step 4: run the init script (depending on windows or unix)