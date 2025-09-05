function [rf, gz, gz_reph] = SIGPY_SLR(flip_angle, pulse_duration, phase_offset, tbw, ptype, ftype, d1, d2, cancel_alpha_phs, dz, sys)

%% calculate pulse object using: sigpy.mri.rf.slr.dzrf
% sigpy-parameters:
%    n       number of time points
%    tbw     pulse time bandwidth product
%    ptype   pulse type:  st  (small-tip excitation)
%                         ex  (pi/2 excitation pulse)
%                         se  (spin-echo pulse)
%                         inv  (inversion)
%                         sat  (pi/2 saturation pulse)
%    ftype   type of filter to use:  ms  (sinc)
%                                    pm  (Parks-McClellan equal-ripple)
%                                    min  (minphase using factored pm)
%                                    max  (maxphase using factored pm)
%                                    ls  (least squares)
%    d1      passband ripple level
%    d2      stopband ripple level
%    cancel_alpha_phs - bool
%    for  ex  pulses, absorb the alpha phase profile from beta s profile,
%    so they cancel for a flatter total phase

n = round(pulse_duration/sys.rfRasterTime);

%% search for sigpy pulse files
sigpy_id = [ 'sigpy_slr_' ...
             ptype '_' ...
             num2str(pulse_duration*1e3, '%.2f') 'ms_' ...
             num2str(n) 'points_' ...
             num2str(tbw, '%.1f') 'tbw_' ...
             ftype '_' ...             
             num2str(d1, '%.3f') 'd1_' ...
             num2str(d2, '%.3f') 'd2_' ...
             num2str(cancel_alpha_phs, '%.0f') 'cap'];

sigpy_id    = [strrep(sigpy_id, '.', ',') '.mat'];
sigpy_path  = pulseq_get_path('SIGPY_find_waveforms');
sigpy_path  = [sigpy_path sigpy_id];
sigpy_exist = isfile(sigpy_path);

if sigpy_exist==1
    load(sigpy_path);
end

%% calculate sigpy pulse shape
if sigpy_exist==0
   
    % string: cancel_alpha_phs
    if cancel_alpha_phs==0
        cap_bool = 'False';
    else
        cap_bool = 'True';
    end

    % read sigpy environment and python specs 
    if ispc()
        temp_py_file = which('init_SIGPY_windows_unix');
        temp_py_file = [temp_py_file(1:end-25) 'py_cmd_windows.txt'];
        temp_py_file = fopen(temp_py_file);
        temp_cmd     = fgetl(temp_py_file);
        fclose(temp_py_file);
    elseif isunix()
        temp_py_file = which('init_SIGPY_windows_unix');
        temp_py_file = [temp_py_file(1:end-25) 'py_cmd_unix.txt'];
        temp_py_file = fopen(temp_py_file);
        temp_cmd     = fgetl(temp_py_file);
        fclose(temp_py_file);
    else
        error('I hate fruits...')
    end

    % command for SLR pulse
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

    % execute command in terminal
    [status, result] = system(cmd);
    if status~=0
        error('executing python command failed');
    end
    
    % read pulse shape: complex waveform [Hz]
    lines  = regexp(result,'\n','split'); 
    signal = str2num(lines{1});
    signal = signal / max(abs(signal));

    % delete peaks
    if strcmp(ftype, 'pm') || strcmp(ftype, 'min') || strcmp(ftype, 'max')
        signal(1) = signal(2);
        signal(n) = signal(n-1);
    end

    % save sigpy file
    save(sigpy_path, 'signal');

end

%% create sigpy pulse object and gradients

% scale waveform for flip angle
signal = signal / (abs(sum(signal)) * sys.rfRasterTime *2*pi) *flip_angle;

% create dummy pulse object
[rf, gz, gz_reph] = mr.makeSincPulse( pi/100, ... % only a dummy
                                      sys, ...
                                      'Duration',       pulse_duration,...
                                      'SliceThickness', dz, ...
                                      'timeBwProduct',  tbw, ...
                                      'PhaseOffset',    0, ...
                                      'use', 'excitation' );

% replace pulse shape and phase offset
rf.signal      = rf.signal*0;
rf.signal(1:n) = signal(:);
rf.phaseOffset = phase_offset;

% rf use
if strcmp(ptype, 'st')
    rf.use = 'excitation';
end
if strcmp(ptype, 'ex')
    rf.use = 'excitation';
end
if strcmp(ptype, 'se')
    rf.use = 'refocusing';
end
if strcmp(ptype, 'inv')
    rf.use = 'inversion';
end
if strcmp(ptype, 'sat')
    rf.use = 'saturation';
end    

end