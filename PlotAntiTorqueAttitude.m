function PlotAntiTorqueAttitude(params, DN, SN_anti, N_L2, N_L2_anti)
% Plot the anti-torque attitude geometry and initial/anti torque vectors.
% This mirrors the inline anti-torque attitude plot in
% Electrostatic_Momentum_Management.m for side-by-side comparison.

plotting_servicer = params.servicer.S_spheres;
for i = 1:size(params.servicer.S_spheres,2)
    plotting_servicer(1:3,i) = SN_anti'*params.servicer.S_spheres(1:3,i);
end

serv_N_COM = SN_anti'*params.servicer.S_COM;
deb_N_COM = DN'*params.debris.D_COM;

if norm(N_L2) > 0
    N_L2_norm = N_L2/norm(N_L2);
else
    N_L2_norm = [0;0;0];
end

if norm(N_L2_anti) > 0
    N_L2_anti_norm = N_L2_anti/norm(N_L2_anti);
else
    N_L2_anti_norm = [0;0;0];
end

figure
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.D_spheres, plotting_servicer,...
    [0 0 0], params.N_rvec_km*1000, [params.debris.voltage, params.servicer.voltage])

p(1)=quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
    serv_N_COM(2) + params.N_rvec_km(2)*1000,...
    serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_norm(1),N_L2_norm(2),N_L2_norm(3),6,'Linewidth',3);
p(2)=quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000,...
    serv_N_COM(2) + params.N_rvec_km(2)*1000,...
    serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_anti_norm(1),N_L2_anti_norm(2),N_L2_anti_norm(3),6,'Linewidth',3);

scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')
colormap(parula)
c=colorbar;
c.Label.String = 'Surface Charge Density (nC/m^2)';
c.Label.FontSize = 14;
xlim([-3,40])
ylim([-17,17])
zlim([-12,12])
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid off
view(3)
hold off
end
