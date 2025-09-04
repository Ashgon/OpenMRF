function [k] = SPI_traj_loop_in(k, t_spoil, tukey_r, system)

    % calculate spoiling trajectory: rewinders
    n_spoil  = round(t_spoil / system.gradRasterTime);
    
    kr       = abs(k);
    kr_spoil = kr(end) + (1:n_spoil) * (kr(end)-kr(end-1));
    my_tukey = 1-tukeywin(2*n_spoil, tukey_r)';
    my_tukey(n_spoil+1:end)=[];
    kr_spoil = kr_spoil .* my_tukey;
    kr_spoil(end) = 0;
    
    kt       = unwrap(atan2(imag(k),real(k)));
    kt_slope = kt(end) - kt(end-1);
    kt_slope = linspace(kt_slope, 0, n_spoil );
    kt_spoil = zeros(1, n_spoil);
    kt_inc   = kt(end);
    for j=1:n_spoil
        kt_inc      = kt_inc + kt_slope(j);
        kt_spoil(j) = kt_inc;
    end
    
    kr = [kr kr_spoil];
    kt = [kt kt_spoil];
    k  = kr .* exp(1i .* kt);

end

