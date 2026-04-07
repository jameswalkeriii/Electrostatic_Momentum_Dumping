function ShowPlots(params,results,ttotal,x0_H,flag) 
% Set plot parameters
set(0,'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.5,'DefaultAxesTitleFontWeight', 'bold') ;

    if flag == true
        time_hours = results.Ttot/3600;
        mode_hist = results.flags(2,:);
        mode_colors = [0.85 0.92 1.00;
                       1.00 0.93 0.80;
                       0.86 0.95 0.86;
                       0.96 0.86 0.86];
                   
        normH = zeros(1,length(results.Htot_servicer));
        for i = 1:length(results.Htot_servicer)
            normH(i) = norm(results.Htot_servicer(:,i));
        end
        
        C2_intial = MRP2C(results.Xtot_servicer(1:3,1));
        
        intial_sph_loc = params.servicer.N_spheres*0;
        for i = 1:length(params.servicer.N_spheres)
            intial_sph_loc(1:3,i) = C2_intial*params.servicer.N_spheres(1:3,i);
            intial_sph_loc(4,i) = params.servicer.N_spheres(4,i);
        end
        C2_end = MRP2C(results.Xtot_servicer(1:3,end));
        final_sph_loc = params.servicer.N_spheres*0;
        for i = 1:length(params.servicer.N_spheres)
            final_sph_loc(1:3,i) = C2_end*params.servicer.N_spheres(1:3,i);
            final_sph_loc(4,i) = params.servicer.N_spheres(4,i);
        end

        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.Xtot_servicer(1:3,:), false))
        plot(time_hours,results.Xtot_servicer(1:3,:),'Linewidth',2)
        xlabel('Time (hours)')
        ylabel('MRP Components')
        xlim([time_hours(1), time_hours(end)])

        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.normBN_servicer, false))
        plot(time_hours,results.normBN_servicer,'Linewidth',2)
        xlabel('Time (hours)')
        ylabel('MRP Magnitude')
        xlim([time_hours(1), time_hours(end)])

        

        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.Xtot_servicer(7:9,:)*60*2*pi, false))
        plot(time_hours,results.Xtot_servicer(7:9,:)*60*2*pi,'Linewidth',2)
        xlabel('Time (hours)','Fontsize',14)
        ylabel('Reaction Wheel Speeds (rpm)','Fontsize',14)
        legend('Om\_1','Om\_2','Om\_3','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])
        
        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.reftot_servicer(1:3,:), false))
        plot(time_hours,results.reftot_servicer(1:3,:),'Linewidth',2)
        xlabel('Time (hours)')
        ylabel('Desired Attitude (MRPs)')
        xlim([time_hours(1), time_hours(end)])

        figure
        subplot(2,1,1)
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.aterrtot_servicer, true))
        semilogy(time_hours,results.aterrtot_servicer,'Linewidth',2)
        ylabel('Magnitude sigma B/R ','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])
        subplot(2,1,2)
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.werrtot_servicer, true))
        semilogy(time_hours,results.werrtot_servicer,'Linewidth',2)
        xlabel('Time (hours)','Fontsize',14)
        ylabel('Magnitude omega B/R ','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])

        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.Xtot_servicer(4:6,:), false))
        plot(time_hours,results.Xtot_servicer(4:6,:),'Linewidth',2)
        xlabel('Time (hours)','Fontsize',14)
        ylabel('Angular Velocity','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])
        
        figure
        mode_change_idx = [1, find(diff(mode_hist) ~= 0) + 1, length(mode_hist) + 1];
        hold on
        for i_mode = 1:length(mode_change_idx)-1
            idx_start = mode_change_idx(i_mode);
            idx_end = mode_change_idx(i_mode+1) - 1;
            active_mode = mode_hist(idx_start);
            x_start = time_hours(idx_start);
            x_end = time_hours(idx_end);
            fill([x_start x_end x_end x_start], [0.5 0.5 4.5 4.5], ...
                mode_colors(active_mode,:), 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        end
        stairs(time_hours,mode_hist,'k','Linewidth',2)
        hold on
        mode_switch_idx = find(diff(mode_hist) ~= 0) + 1;
        for i_mode = 1:length(mode_switch_idx)
            xline(time_hours(mode_switch_idx(i_mode)),'--')
        end
        ylim([0.5,4.5])
        yticks(1:4)
        yticklabels({'First Attitude','Slew to Second','Second Attitude','Slew to First'})
        xlabel('Time (hours)','Fontsize',14)
        ylabel('Desaturation Mode','Fontsize',14)
        title('Desaturation Mode History','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])
        
        figure 
        hold on
        addModeBands(gca, results.Ttot/3600/24, mode_hist, mode_colors, calcPlotLimits(results.Xtot_servicer(7:9,:)*60/(2*pi), false))
        plot(results.Ttot/3600/24,results.Xtot_servicer(7,:)*60/(2*pi),'Linewidth',2)
        plot(results.Ttot/3600/24,results.Xtot_servicer(8,:)*60/(2*pi),'Linewidth',2)
        plot(results.Ttot/3600/24,results.Xtot_servicer(9,:)*60/(2*pi),'color',plot_colors.col_4,'Linewidth',2)
        xlabel('Time (days)','Fontsize',14,'Fontname','Times New Roman')
        ylabel('RW Velocity (rpm)','Fontsize',14,'Fontname','Times New Roman')
        legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','Fontsize',12,'Fontname','Times New Roman','Location','southeast')
        xlim([0,results.Ttot(end)/3600/24])
        
        
        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(vecnorm(results.Htot_servicer,2,1), false))
        plot(time_hours,vecnorm(results.Htot_servicer,2,1),'Linewidth',2)
        ylabel('Angular Momentum (Nms)','Fontsize',14,'Fontname','Times New Roman')
        xlabel('Time (hours)','Fontsize',14,'Fontname','Times New Roman')
        xlim([time_hours(1), time_hours(end)])

        figure
        hold on
        addModeBands(gca, time_hours, mode_hist, mode_colors, calcPlotLimits(results.L_e_servicer_tot, false))
        plot(time_hours,results.L_e_servicer_tot,'Linewidth',2)
        ylabel('Electrostatic Torques (Nm)','Fontsize',14)
        xlabel('Time (hours)','Fontsize',14)
        legend('L\_x','L\_y','L\_z','Fontsize',14)
        xlim([time_hours(1), time_hours(end)])

        %Initial position plot
        figure 
        grid on
        hold on
        xlabel('Relative X distance [m]')
        ylabel('Relative Y distance [m]')
        zlabel('Relative Z distance [m]')
        set(gca,'FontName','times')
        makeSphsPicture_2craft( params.debris.N_spheres, intial_sph_loc, x0_H*1000,...
        [0 0 0], params.V)
%         quiver3(0,0,0,4*Htot(1,1)/100,4*Htot(2,1)/100,4*Htot(3,1)/100,'r','Linewidth',3)
%         quiver3(0,0,0,10e3*L_e_tot(1,1),10e3*L_e_tot(2,1),10e3*L_e_tot(3,1),'b','Linewidth',3)
%         quiver3(0,0,0,10e3*L_e_tot(1,11051),10e3*L_e_tot(2,11051),10e3*L_e_tot(3,11051),'b','Linewidth',3)
        axis tight
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
%         c.Label.FontSize = 14;
%         title('Initial Positions')
        view(3)
        
        %Final position plot
        figure 
        grid on
        hold on
        xlabel('Relative X distance [m]')
        ylabel('Relative Y distance [m]')
        zlabel('Relative Z distance [m]')
        set(gca,'FontName','times')
        makeSphsPicture_2craft( params.debris.N_spheres, final_sph_loc, x0_H*1000,...
        [0 0 0], params.V)
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
%         c.Label.FontSize = 14;
        axis tight
%         title('Final Positions')
        view(3)
    end

end

function addModeBands(ax, t_hist, mode_hist, mode_colors, y_limits)
    mode_change_idx = [1, find(diff(mode_hist) ~= 0) + 1, length(mode_hist) + 1];
    for i_mode = 1:length(mode_change_idx)-1
        idx_start = mode_change_idx(i_mode);
        idx_end = mode_change_idx(i_mode+1) - 1;
        active_mode = mode_hist(idx_start);
        x_start = t_hist(idx_start);
        x_end = t_hist(idx_end);
        patch(ax, [x_start x_end x_end x_start], [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
            mode_colors(active_mode,:), 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    ylim(ax, y_limits)
    set(ax, 'Layer', 'top')
end

function y_limits = calcPlotLimits(data, use_log)
    vals = data(:);
    vals = vals(isfinite(vals));
    if use_log
        vals = vals(vals > 0);
    end

    if isempty(vals)
        y_limits = [0.1 1];
        return
    end

    y_min = min(vals);
    y_max = max(vals);

    if y_min == y_max
        if use_log
            y_limits = [y_min/10, y_max*10];
        else
            pad = max(1, abs(y_min)*0.1 + 1e-6);
            y_limits = [y_min-pad, y_max+pad];
        end
        return
    end

    if use_log
        y_limits = [10^floor(log10(y_min)), 10^ceil(log10(y_max))];
    else
        pad = 0.05*(y_max-y_min);
        y_limits = [y_min-pad, y_max+pad];
    end
end
