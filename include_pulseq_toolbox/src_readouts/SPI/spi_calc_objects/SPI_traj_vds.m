function [traj] = SPI_traj_vds(SPI, FOV, system, mode_plot)

    if nargin<4
        mode_plot = 0;
    end
    
    %% SI -> hargreaves-units
    Nvds     = SPI.geo.Nvds;
    t_dwell  = SPI.geo.t_dwell;
    t_dwell  = round( t_dwell / system.adcRasterTime ) * system.adcRasterTime;
    Fcoeff   = SPI.geo.Fcoeff * FOV.fov_xy *1e2; % [cm]
    rmax     = SPI.geo.kmax *1e-2; % [1/cm]              
    grad_lim = SPI.geo.grad_lim;
    slew_lim = SPI.geo.slew_lim;
    max_grad = grad_lim * system.maxGrad/system.gamma *1e2; % [T/cm]
    max_slew = slew_lim * system.maxSlew/system.gamma *1e2; % [T/cm/s]
    [k, g, s, t, ~, ~] = vds(max_slew, max_grad, t_dwell, Nvds, Fcoeff, rmax);
    
    %% SI
    k = k *1e2; % SI
    g = g *1e-2; % SI
    s = s *1e-2; % SI
    
    %% disp params
    disp(' ');
    disp('----------------------- TRAJ ----------------------- ');
    disp(['N samples per interleave: ' num2str(numel(t)) ' (' num2str(round(1/t_dwell/1000)) 'kHz)']);
    disp(['adc time: ' num2str(t(end)*1000, '%.1f') 'ms']);
    disp(['Undersampling for single interleave: '  num2str(numel(t) / (FOV.Nxy^2*pi/4) *100, '%.1f') '%' ]);
    disp(['Undersampling for all interleaves: '  num2str(Nvds * numel(t) / (FOV.Nxy^2*pi/4) *100, '%.1f') '%' ]);
    disp('--------------------------------------------------- ');
    disp(' ');
    
    %% plot: kspace, gradients, slew rates, dcf2d
    if mode_plot==1
    
        % DCF
        temp_k   = zeros(Nvds, numel(k));
        temp_phi = linspace(0, 2*pi, Nvds+1);
        temp_phi(end) = [];
        for j=1:Nvds
            temp_k(j,:) = k * exp(1i * temp_phi(j));
        end
        temp_kx    = real(temp_k(:));
        temp_ky    = imag(temp_k(:));
        temp_N     = [1 1] * FOV.Nxy;
        temp_fov   = FOV.fov_xy;
        temp_nufft = {temp_N, [6 6], 2*temp_N, temp_N/2, 'table', 2^12, 'minmax:kb'};
        temp_G     = Gmri([temp_kx, temp_ky], true(temp_N), 'fov', temp_fov, 'basis', {'rect'}, 'nufft', temp_nufft);
        temp_dcf   = abs(mri_density_comp([temp_kx, temp_ky], 'pipe', 'G', temp_G.arg.Gnufft));
        DCF2D      = mg_get_heatmap(temp_kx, temp_ky, temp_dcf, 2);
    
        figure()
        subplot(2,3,1);
        plot(real(k),imag(k), 'k.');
        title('ky vs kx');
        axis('square');
     
        subplot(2,3,2);
        plot(real(temp_k(:)),imag(temp_k(:)), 'k.');
        title('ky vs kx: full');
        axis('square');
    
        subplot(2,3,3)
        imagesc(DCF2D); axis image; axis off; colormap(plasma(1000));
        title('DCF2d');
    
        subplot(2,3,4);
        plot(t*1e3,real(k),'r-', t*1e3,imag(k),'b-', t*1e3,abs(k),'k-');
        title('k vs t');
        xlabel('time [ms]')
        ylabel('k (m^{-1})');
        
        subplot(2,3,5);
        plot(t*1e3,real(g)*1e3,'r-', t*1e3,imag(g)*1e3,'b-', t*1e3,abs(g)*1e3,'k-');
        title('g vs t');
        xlabel('time [ms]')
        ylabel('G (mT/m)');
        
        subplot(2,3,6);
        plot(t*1e3,real(s),'r-', t*1e3,imag(s),'b-', t,abs(s),'k-');
        title('s vs t');
        xlabel('time [ms]')
        ylabel('Slew Rate (T/m/s)');
    
    end
    
    %% interpolate vds trajectory to grad raster
    tq     = (0 : 1 : ceil(t(end)/system.gradRasterTime)) *  system.gradRasterTime;
    k      = interp1(t, k, tq);
    k(1)   = 0;
    k(end) = [];
    
    %% output vds trajectory
    traj.mode      = 'vds';
    traj.k         = k;
    traj.n_samples = numel(t);
    traj.t_adc     = t(end);
    traj.t_dwell   = t_dwell;
    traj.bandwidth = 1/traj.t_dwell;

end

