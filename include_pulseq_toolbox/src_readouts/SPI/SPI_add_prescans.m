% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% -------------------- Noise pre-scans --------------------
% ---------------------------------------------------------

if isfield(SPI, 'Nnoise')
    for loop_noise = 1:SPI.Nnoise
        seq.addBlock(SPI.adc, ceil((mr.calcDuration(SPI.adc) + 1e-3) / system.blockDurationRaster) * system.blockDurationRaster);
    end
end
