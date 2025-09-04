%% add sequence objects for cardio or abdominal MRF

% init loop counters for contrast preparations
loop_SAT    = 1;
loop_INV    = 1;
loop_T2     = 1;
loop_SL     = 1;
loop_MLEV   = 1;
loop_ADIASL = 1;

% reset rf spoiling increment
loop_rf_inc = 0;  
   
for loop_MRF = 1 : MRF.n_segm

	% Segment label for GE scanners
	if flag_GE==1
		seq.addBlock(mr.makeDelay(system.gradRasterTime), mr.makeLabel('SET', 'TRID', loop_MRF+1));
	end
 
    % Trigger: R-Wave
    if exist('TRIG_IN', 'var')
        seq.addBlock(TRIG_IN);   
    end

    % Soft Delay and Dynamic Delay after Trigger (adjust cardiac phase)
    if strcmp(MRF.delay_soft.type, 'softDelay')
        seq.addBlock(system.gradRasterTime, MRF.delay_soft);
        seq.addBlock(MRF.delay_dynamic(loop_MRF));
    else
        seq.addBlock(mr.makeDelay(mr.calcDuration(MRF.delay_soft) + mr.calcDuration(MRF.delay_dynamic(loop_MRF))));
    end    

    % MRF Preparations:
    MRF_add_preparation();

    % Spiral Readouts
    for loop_seg = 1:SPI.nr
        loop_NR = (loop_MRF-1)*SPI.nr + loop_seg;
        SPI_add();
    end

end