function ShowPlots(Xtot, Htot, Ttot, ttot, aterrtot, werrtot,...
    ustot,Lrtot,L_e_tot,params0,params,normBN,reftot,ttotal,x0_H,flag) 
    col = [.466,.674,.188];
    if flag == true
        
        for i = 1:length(Htot)
            normH(i) = norm(Htot(:,i));
        end
        
        C2_intial = MRP2C(Xtot(1:3,1));
        
        intial_sph_loc = params.servicer.spheres*0;
        for i = 1:length(params.servicer.spheres)
            intial_sph_loc(1:3,i) = C2_intial*params.servicer.spheres(1:3,i);
            intial_sph_loc(4,i) = params.servicer.spheres(4,i);
        end
        C2_end = MRP2C(Xtot(1:3,end));
%         C2_end = MRP2C([-1;0;0]);
        final_sph_loc = params.servicer.spheres*0;
        for i = 1:length(params.servicer.spheres)
            final_sph_loc(1:3,i) = C2_end*params.servicer.spheres(1:3,i);
            final_sph_loc(4,i) = params.servicer.spheres(4,i);
        end

        figure
        plot(ttot,Xtot(1:3,:),'Linewidth',2)
        xlabel('Time (s)')
        ylabel('MRP Components')

        figure
        plot(ttot(2:end),normBN,'Linewidth',2)
        xlabel('Time (s)')
        ylabel('MRP Magnitude')

        

        figure
        plot(ttot,Xtot(7:9,:)*60*2*pi,'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Reaction Wheel Speeds (rpm)','Fontsize',14)
        legend('Om\_1','Om\_2','Om\_3','Fontsize',14)
        xlim([0,ttotal])
        
        figure
        plot(ttot(2:end),reftot(1:3,:),'Linewidth',2)
        xlabel('Time (s)')
        ylabel('Desired Attitude (MRPs)')

        figure
        subplot(2,1,1)
        semilogy(ttot(2:end),aterrtot,'Linewidth',2)
        ylabel('Magnitude sigma B/R ','Fontsize',14)
        subplot(2,1,2)
        semilogy(ttot(2:end),werrtot,'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Magnitude omega B/R ','Fontsize',14)
        
        figure 
        plot(ttot/3600/24,Xtot(7,:)*60*2*pi,'Linewidth',2)
        hold on
        plot(ttot/3600/24,Xtot(8,:)*60*2*pi,'Linewidth',2)
        plot(ttot/3600/24,Xtot(9,:)*60*2*pi,'color',col,'Linewidth',2)
        xlabel('Time (days)','Fontsize',14,'Fontname','Times New Roman')
        ylabel('RW Velocity (rpm)','Fontsize',14,'Fontname','Times New Roman')
        legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','Fontsize',12,'Fontname','Times New Roman','Location','southeast')
        xlim([0,ttotal/3600/24])
        
        figure
        plot(ttot(2:end),normH,'Linewidth',2)
        hold on
        xline(38120,'--')
        ylabel('Angular Momentum (Nms)','Fontsize',14,'Fontname','Times New Roman')
        xlabel('Time (s)','Fontsize',14,'Fontname','Times New Roman')

        figure
        plot(ttot(2:end),L_e_tot,'Linewidth',2)
        ylabel('Electrostatic Torques (Nm)','Fontsize',14)
        xlabel('Time (s)','Fontsize',14)
        legend('L\_x','L\_y','L\_z','Fontsize',14)

        %Initial position plot
        figure 
        grid on
        hold on
        xlabel('Relative X distance [m]')
        ylabel('Relative Y distance [m]')
        zlabel('Relative Z distance [m]')
        set(gca,'FontName','times')
        makeSphsPicture_2craft( params.debris.spheres, intial_sph_loc, x0_H*1000,...
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
        makeSphsPicture_2craft( params.debris.spheres, final_sph_loc, x0_H*1000,...
        [0 0 0], params.V)
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
%         c.Label.FontSize = 14;
        axis tight
%         title('Final Positions')
        view(3)
    end

end