function PlotInitialPosition(params, SN, DN, N_L2)
% Plot the initial target/servicer geometry and the initial servicer torque.
% This mirrors the inline initial-position plot in
% Electrostatic_Momentum_Management.m for side-by-side comparison.

plotting_servicer = params.servicer.S_spheres;
for i = 1:size(params.servicer.S_spheres,2)
    plotting_servicer(1:3,i) = SN'*params.servicer.S_spheres(1:3,i);
end

plotting_debris = params.debris.D_spheres;
for i = 1:size(params.debris.D_spheres,2)
    plotting_debris(1:3,i) = DN'*params.debris.D_spheres(1:3,i);
end

serv_N_COM = SN'*params.servicer.S_COM;
deb_N_COM = DN'*params.debris.D_COM;

if norm(N_L2) > 0
    N_L2_norm = N_L2/norm(N_L2);
else
    N_L2_norm = [0;0;0];
end

figure
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(plotting_debris, plotting_servicer,...
    [0 0 0], params.N_rvec_km*1000,  [params.debris.voltage, params.servicer.voltage])

quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
    serv_N_COM(2) + params.N_rvec_km(2)*1000,...
    serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_norm(1),N_L2_norm(2),N_L2_norm(3),6,'Linewidth',3)

scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')

% 
% c=colorbar;
% c.Label.String = 'Surface Charge Density (nC/m^2)';
% c.Label.FontSize = 14;

axis equal
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
