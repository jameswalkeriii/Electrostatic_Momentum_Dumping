function PlotServicerInitialAntiOverlayContinuous(params, DN, SN_initial, SN_anti, N_L_initial, N_L_anti)
% Plot target plus initial and anti-torque servicer attitudes as continuous
% surfaces. Each electrostatic torque is drawn from its respective servicer
% center of mass with a common scale factor.

N_servicer_origin = params.N_rvec_km*1000;
N_initial_COM = N_servicer_origin + SN_initial'*params.servicer.S_COM;
N_anti_COM = N_servicer_origin + SN_anti'*params.servicer.S_COM;
N_debris_COM = DN'*params.debris.D_COM;

max_torque = max([norm(N_L_initial), norm(N_L_anti)]);
if max_torque == 0
    max_torque = 1;
end
torque_scale = 7/max_torque;

figure
hold on
set(gca,'FontName','times')

plotTargetSurfaceLocal([0;0;0], DN, [0.72 0.76 0.80], 0.50)
plotServicerSurfaceLocal(N_servicer_origin, SN_initial, [0.10 0.32 0.75], 0.34)
plotServicerSurfaceLocal(N_servicer_origin, SN_anti, [0.72 0.14 0.48], 0.34)

p(1) = plotTorqueLocal(N_initial_COM, N_L_initial, torque_scale, [0.05 0.22 0.85], 'Initial attitude torque');
p(2) = plotTorqueLocal(N_anti_COM, N_L_anti, torque_scale, [0.78 0.08 0.18], 'Anti-torque attitude torque');

scatter3(N_debris_COM(1),N_debris_COM(2),N_debris_COM(3),20,'r','filled','HandleVisibility','off')
scatter3(N_initial_COM(1),N_initial_COM(2),N_initial_COM(3),20,'k','filled','HandleVisibility','off')
scatter3(N_anti_COM(1),N_anti_COM(2),N_anti_COM(3),20,'k','filled','HandleVisibility','off')

axis equal
xlim([-3,40])
ylim([-17,17])
zlim([-12,12])
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
% legend(p,'Location','best')
grid off
view(3)
camlight headlight
lighting gouraud
hold off
end

function h = plotTorqueLocal(origin, torque, torque_scale, color, display_name)
if norm(torque) == 0
    h = plot3(nan,nan,nan,'Color',color,'DisplayName',display_name);
    return
end
scaled_torque = torque_scale*torque;
h = quiver3(origin(1), origin(2), origin(3),...
    scaled_torque(1), scaled_torque(2), scaled_torque(3), 0,...
    'LineWidth',3,'Color',color,'MaxHeadSize',0.65,'DisplayName',display_name);
end

function plotServicerSurfaceLocal(offset, SN, face_color, face_alpha)
bus_limits = [0.55, 5.55;
              -2.5, 2.5;
              -2.5, 2.5];

panel_limits = [1.05, 1.95;
                -12.6, -2.5;
                -2.4, 2.4];

drawBoxSurfaceLocal(bus_limits, offset, SN, face_color, face_alpha)
drawBoxSurfaceLocal(panel_limits, offset, SN, face_color, face_alpha)
end

function plotTargetSurfaceLocal(offset, DN, face_color, face_alpha)
[bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal();
drawBoxSurfaceLocal(bus_limits, offset, DN, face_color, face_alpha)
drawBoxSurfaceLocal(panel_limits_1, offset, DN, face_color, face_alpha)
drawBoxSurfaceLocal(panel_limits_2, offset, DN, face_color, face_alpha)
end

function [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal()
% Hard-coded approximate SSL-1300 target geometry in the same rotated frame
% as params.debris.D_spheres.
bus_limits = [0.6, 3.4;
              -1.05, 1.05;
              -1.0, 1.0];

panel_limits_1 = [0.625, 1.375;
                  1.05, 15.7;
                  -1.3, 1.8];

panel_limits_2 = [0.625, 1.375;
                  -15.7, -1.05;
                  -1.3, 1.8];
end

function drawBoxSurfaceLocal(limits, offset, BN, face_color, face_alpha)
x_min = limits(1,1);
x_max = limits(1,2);
y_min = limits(2,1);
y_max = limits(2,2);
z_min = limits(3,1);
z_max = limits(3,2);

vertices_B = [x_min y_min z_min;
              x_max y_min z_min;
              x_max y_max z_min;
              x_min y_max z_min;
              x_min y_min z_max;
              x_max y_min z_max;
              x_max y_max z_max;
              x_min y_max z_max];
vertices_N = (BN'*vertices_B')' + offset(:)';
faces = [1 2 3 4;
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];
patch('Vertices',vertices_N,'Faces',faces,...
    'FaceColor',face_color,'FaceAlpha',face_alpha,...
    'EdgeColor',[0.08 0.10 0.12],'LineWidth',0.8,...
    'HandleVisibility','off')
end
