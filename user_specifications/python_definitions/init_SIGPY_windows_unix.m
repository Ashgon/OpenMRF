%% init on windows
clear

% try simple SLR pulse
n = 1000;
tbw = 6;
ptype = 'st';
ftype = 'ls';
d1 = 0.01;
d2 = 0.01;
cap_bool = 'False';

% enter your sigpy environment and python specs
temp_env = 'C:\\ProgramData\\miniforge3\\Scripts\\activate.bat SIGPY';
temp_py  = 'C:\\Users\\mag28ev\\.conda\\envs\\SIGPY\\python';

% save in .txt file
temp_cmd = [ 'cmd /c ' temp_env ' && ' temp_py ' ']; 
writelines(temp_cmd, 'py_cmd_windows.txt');

% create command
cmd = [ temp_cmd ...
        '-c "import sigpy.mri.rf; ' ...
        'mypulse = sigpy.mri.rf.dzrf( ' ...
        num2str(n) ', ' ... 
        num2str(tbw) ', ' ...
        'ptype = ''' ptype ''', ' ...
        'ftype = ''' ftype ''', ' ...
        'd1 = ' num2str(d1) ', ' ...
        'd2 = ' num2str(d2) ', ' ...
        'cancel_alpha_phs = ' cap_bool ...
        ' ); print(*mypulse)"' ];

% execute command and read pulse waveform
[status, result] = system(cmd);
lines  = regexp(result,'\n','split'); 
signal = str2num(lines{1});
signal = signal / max(abs(signal));

figure()
plot(abs(signal))

%% init on unix
clear

% try simple SLR pulse
n = 1000;
tbw = 6;
ptype = 'st';
ftype = 'ls';
d1 = 0.01;
d2 = 0.01;
cap_bool = 'False';

% enter your sigpy environment and python specs
temp_env = '/home/ubuntu/miniforge3/bin/activate SIGPY';
temp_py  = 'python3';

% save in .txt file
temp_cmd = [ 'source ' temp_env ' && ' temp_py ' ']; 
writelines(temp_cmd, 'py_cmd_unix.txt');

% create command
cmd = [ temp_cmd ...
        '-c "import sigpy.mri.rf; ' ...
        'mypulse = sigpy.mri.rf.dzrf( ' ...
        num2str(n) ', ' ... 
        num2str(tbw) ', ' ...
        'ptype = ''' ptype ''', ' ...
        'ftype = ''' ftype ''', ' ...
        'd1 = ' num2str(d1) ', ' ...
        'd2 = ' num2str(d2) ', ' ...
        'cancel_alpha_phs = ' cap_bool ...
        ' ); print(*mypulse)"' ];


% execute command and read pulse waveform
[status, result] = system(cmd);
lines  = regexp(result,'\n','split'); 
signal = str2num(lines{1});
signal = signal / max(abs(signal));

figure()
plot(abs(signal))