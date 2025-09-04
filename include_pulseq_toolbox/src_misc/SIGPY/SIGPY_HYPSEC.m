function [rf] = SIGPY_HYPSEC(beta, mu, pulse_duration, phase_offset, sys)

%% calculate pulse object using: sigpy.mri.rf.adiabatic.hypsec
% sigpy-parameters:
%    n       number of time points
%    beta  - AM waveform parameter
%    mu    - a constant, determines amplitude of frequency sweep

n = round(pulse_duration/sys.rfRasterTime);

%% search for sigpy pulse files
sigpy_id = [ 'sigpy_hypsec_' ...
             num2str(pulse_duration*1e3, '%.2f') 'ms_' ...
             num2str(n) 'points_' ...
             num2str(beta, '%.2f') 'beta_' ...
             num2str(mu, '%.2f') 'mu'];

sigpy_id    = [strrep(sigpy_id, '.', ',') '.mat'];
sigpy_path  = pulseq_get_path('SIGPY_find_waveforms');
sigpy_path  = [sigpy_path sigpy_id];
sigpy_exist = isfile(sigpy_path);

if sigpy_exist==1
    load(sigpy_path);
end

%% calculate sigpy pulse shape

if sigpy_exist==0
    
    % read sigpy environment and python specs 
    if ispc()
        temp_py_file = which('SIGPY_HYPSEC');
        temp_py_file = [temp_py_file(1:end-14) 'py_cmd_windows.txt'];
        temp_py_file = fopen(temp_py_file);
        temp_cmd     = fgetl(temp_py_file);
        fclose(temp_py_file);
    elseif isunix()
        temp_py_file = which('SIGPY_HYPSEC');
        temp_py_file = [temp_py_file(1:end-14) 'py_cmd_unix.txt'];
        temp_py_file = fopen(temp_py_file);
        temp_cmd     = fgetl(temp_py_file);
        fclose(temp_py_file);
    else
        error('I hate fruits...')
    end

    % commands for HYPSEC pulse
    cmd1 = [ temp_cmd ...
             '-c "import sigpy.mri.rf; ' ...
             'mypulse = sigpy.mri.rf.adiabatic.hypsec( ' ...
             num2str(n) ', ' ... 
             num2str(beta) ', ' ...
             num2str(mu) ', ' ...
             num2str(pulse_duration) ...
             ' ); print(*mypulse[0]);"' ];
    cmd2 = [ temp_cmd ...
             '-c "import sigpy.mri.rf; ' ...
             'mypulse = sigpy.mri.rf.adiabatic.hypsec( ' ...
             num2str(n) ', ' ... 
             num2str(beta) ', ' ...
             num2str(mu) ', ' ...
             num2str(pulse_duration) ...
             ' ); print(*mypulse[1]);"' ];

    % execute command in terminal
    [status, result1] = system(cmd1);
    if status~=0
        error('executing python command failed');
    end
    [status, result2] = system(cmd2);
    if status~=0
        error('executing python command failed');
    end

    % read pulse shape: ampitude modulation
    lines     = regexp(result1,'\n','split'); 
    signal_am = str2num(lines{1}); % [0...1]
    signal_am = signal_am / max(signal_am) * beta;    
    % read pulse shape: offresonance modulation
    lines     = regexp(result2,'\n','split'); 
    signal_fm = str2num(lines{1}); % [rad/s]
    
    % calculate phase modulation
    signal_phase = cumsum(signal_fm) * sys.rfRasterTime;
    
    % calculate complex pulse shape
    signal = signal_am .* exp(1i* signal_phase);
    
    % save sigpy file
    save(sigpy_path, 'signal');    

end

%% create sigpy pulse object

% create dummy pulse object
rf = mr.makeSincPulse( pi, ...
                       sys, ...
                       'Duration', pulse_duration,...
                       'SliceThickness', 1e6, ...
                       'PhaseOffset', 0, ...
                       'use', 'inversion' );

% replace pulse shape and phase offset
rf.signal      = rf.signal*0;
rf.signal(1:n) = signal(:);
rf.phaseOffset = phase_offset;

end

