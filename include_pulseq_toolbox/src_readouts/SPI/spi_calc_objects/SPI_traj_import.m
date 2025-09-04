function [traj] = SPI_traj_import(SPI, system)

    %% load complex gxy(t) or kxy(t) and duration
    load(SPI.geo.path);
    if exist('gxy', 'var')
        SPI.import.gxy     = gxy;
        SPI.import.gxy_dur = dur;
    elseif exist('kxy', 'var')
        SPI.import.kxy     = kxy;
        SPI.import.kxy_dur = dur;
        gxy                = [0; diff(kxy)];
        gxy                = gxy / max(abs(gxy));
    end
    
    %% interpolate imported trajectory to grad raster time
    t        = linspace(0, dur, numel(gxy))';
    tq       = (0 : 1 : round(dur/system.gradRasterTime))' *  system.gradRasterTime;
    gxy      = interp1(t, gxy, tq);
    k        = cumsum(gxy);
    k(1)     = 0;
    k(end+1) = 0;
    k        = k / max(abs(k));
    k        = k * SPI.geo.kmax;
    
    %% output trajectory
    traj.mode      = 'import';
    traj.k         = k;
    traj.n_samples = ceil(t(end)/SPI.geo.t_dwell);
    traj.t_adc     = t(end);
    traj.t_dwell   = SPI.geo.t_dwell;
    traj.bandwidth = 1/traj.t_dwell;

end

