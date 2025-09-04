function [adc, adcTime, adcNSamples, adcBW, adcDwell] = SPI_calc_adc(SPI, system)

    %% calc adc
    adcTimeDesired     = SPI.geo.traj.t_adc;
    adcNSamplesDesired = SPI.geo.traj.n_samples;
    adcBWDesired       = SPI.geo.traj.bandwidth;
    adcTime            = round( adcTimeDesired / system.blockDurationRaster ) * system.blockDurationRaster - system.blockDurationRaster;
    adcNSamples        = ceil(adcNSamplesDesired/64)*64;
    adcDwell           = round( adcTime / adcNSamples / system.adcRasterTime ) * system.adcRasterTime;
    adcBW              = 1/adcDwell;
    adc                = mr.makeAdc(adcNSamples, 'Dwell', adcDwell); 
    adc.delay          = system.adcDeadTime;
    adc.phaseOffset    = 0; % changed during rf spoiling

    %% prevent adc sampling truncation 
    if adcNSamples > 16384
        error('maximum number of adc samples is 16384! decrease bandwidth!');
    end
    
    %% disp params
    disp(' ');
    disp('----------------------- ADC ----------------------- ');
    disp([ 'N samples desired: '  num2str(adcNSamplesDesired)  '     N samples used: '           num2str(adcNSamples)                ]);
    disp([ 'adc time desired: '   num2str(adcTimeDesired*1e3,  '%.1f') 'ms     adc time used: '  num2str(adcTime*1e3,  '%.1f') 'ms'  ]);
    disp([ 'adc BW desired: '     num2str(adcBWDesired/1e3,    '%.1f') 'kHz    adc BW used: '    num2str(adcBW/1e3,    '%.1f') 'kHz' ]);
    disp('--------------------------------------------------- ');
    disp(' ');

end

