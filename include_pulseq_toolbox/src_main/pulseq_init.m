%% initilization script

% default flags
if ~exist('flag_backup', 'var')
    flag_backup = 0; % 0: off,  1: only backup,  2: backup and send .seq
end
if ~exist('flag_report', 'var')
    flag_report = 0; % 0: off,  1: only timings, 2: full report (slow)
end
if ~exist('flag_pns', 'var')
    flag_pns = 0; % 0: off,  1: simulate PNS stimulation
end
if ~exist('flag_sound', 'var')
    flag_sound = 0; % 0: off,  1: simulate gradient sound
end
if ~exist('flag_mrf', 'var')
    flag_mrf = 0; % 0: off,  1: simulate sequence via MRF toolbox
end
if ~exist('flag_GE', 'var')
    flag_GE = 0; % 0: off, 1: GE segment TRID labels
end

% load user information
if ~exist('pulseq_scanner', 'var')
    pulseq_scanner = [];
end
[pulseq_user, pulseq_lab, pulseq_path, pulseq_scanner, system] = pulseq_get_user_definitions(pulseq_scanner, 1);

% activate GE TRID labels
if strcmp(pulseq_scanner(1:2), 'GE')
    flag_GE = 1;
end

% init seq object
seq = mr.Sequence(system);

% create pulseq ID for linking rawdata with backup of metadata
if ~exist('seq_name', 'var')
    seq_name = 'external';
end
[seq_id, scan_id, wip_id] = pulseq_get_seq_id();  % [yy mm dd tttt] identifier
seq_name = [seq_id '_' pulseq_user '_' seq_name]; % seq_name is used for linking rawdata to backup files