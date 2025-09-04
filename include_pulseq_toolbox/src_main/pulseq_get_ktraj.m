function [ktraj_adc, ktraj_full] = pulseq_get_ktraj(seq, flag_plot)

    warning('OFF', 'mr:restoreShape');
    [ktraj_adc, ~, ktraj_full] = seq.calculateKspacePP();
    if flag_plot==1
        figure();
        hold on
        plot( ktraj_full(1,:), ktraj_full(2,:), 'b-');
        plot( ktraj_adc(1,:),  ktraj_adc(2,:), 'r.');
        axis('equal');
        xlabel('kx [1/m]');
        ylabel('ky [1/m]');
        title('k-space trajectory');
    end

end