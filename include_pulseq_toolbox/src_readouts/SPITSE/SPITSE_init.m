function [SPITSE, ktraj_adc, ktraj_full, ktraj_reco] = SPITSE_init(SPITSE, FOV, system)
    
% ---------------------------------------------------------
% ---------- init parameters and pulseq objects -----------
% ------- readout: SPIral-Turbo-Spin-Echo (SPITSE) --------
% ---------------------------------------------------------

% source: https://github.com/HennigJue/single-shot-spiral-TSE

[seq, SPITSE]           = SPITSE_add([], system, FOV, SPITSE, 1);  % dummy run
[ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, 1);                % calc trajectory
ktraj_reco = ktraj_adc(1:2,:);
SPITSE.plotflag = '000';

end