function PlotServicerTorqueCloud(params, relative_orientation_EMM_Torques)
% Plot all electrostatic torques experienced by the servicer.
% This function mirrors the inline plotting block in
% Electrostatic_Momentum_Management.m so both methods can be compared.

figure
hold on
set(gca,'FontName','times')

plotMSMSphereBodyLocal(params.debris.D_spheres, [0 0 0])
plotServicerBusPanelSurfaceLocal(params.servicer.S_spheres, params.N_rvec_km*1000)

torque_origin = params.servicer.S_COM + params.N_rvec_km*1000;
scatter3(torque_origin(1),torque_origin(2),torque_origin(3),20,'k','filled')

torque_norms = zeros(length(relative_orientation_EMM_Torques),1);
for i_L = 1:length(relative_orientation_EMM_Torques)
    torque_norms(i_L) = norm(relative_orientation_EMM_Torques{i_L}.N_L_elect_serv);
end
max_torque = max(torque_norms);
if max_torque == 0
    max_torque = 1;
end
torque_arrow_scale = 8/max_torque;
torque_colors = torqueMagnitudeColormapLocal(256);

for i_L = 1:length(relative_orientation_EMM_Torques)
    E_torque = relative_orientation_EMM_Torques{i_L}.N_L_elect_serv;
    if torque_norms(i_L) > 0
        color_idx = min(256,max(1,ceil(256*torque_norms(i_L)/max_torque)));
        E_torque = torque_arrow_scale*E_torque;
        quiver3(torque_origin(1), torque_origin(2), torque_origin(3),...
            E_torque(1),E_torque(2),E_torque(3),0,...
            'Linewidth',0.8,'Color',torque_colors(color_idx,:))
    end
end

axis equal
xlim([31,40])
ylim([-17,10])
zlim([-10,10])
colormap(torque_colors)
caxis([0 max_torque])
cb = colorbar;
cb.Label.String = 'Torque Magnitude [N m]';
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid off
view(35,20)
camlight headlight
lighting gouraud
hold off
end

function plotMSMSphereBodyLocal(SPHS, offset)
    n_spheres = size(SPHS,2);
    [xs,ys,zs] = sphere(18);

    for i_sphere = 1:n_spheres
        radius = SPHS(4,i_sphere);
        xp = radius*xs + SPHS(1,i_sphere) + offset(1);
        yp = radius*ys + SPHS(2,i_sphere) + offset(2);
        zp = radius*zs + SPHS(3,i_sphere) + offset(3);
        surf(xp, yp, zp, 'FaceColor', [0.80 0.83 0.86], 'FaceAlpha', 0.35, 'EdgeColor', 'none')
    end
end

function cmap = torqueMagnitudeColormapLocal(n_colors)
    color_stops = [0.10 0.22 0.55;
                   0.00 0.58 0.72;
                   0.18 0.72 0.42;
                   0.74 0.20 0.62;
                   0.45 0.05 0.18];
    stop_locations = linspace(0,1,size(color_stops,1));
    query_locations = linspace(0,1,n_colors);
    cmap = interp1(stop_locations, color_stops, query_locations);
end

function plotServicerBusPanelSurfaceLocal(~, offset)
    bus_limits = [0.55, 5.55;
                  -2.5, 2.5;
                  -2.5, 2.5];

    panel_limits = [1.05, 1.95;
                    -12.6, -2.5;
                    -2.4, 2.4];

    drawBoxSurfaceLocal(bus_limits, offset, [0.72 0.75 0.78], 0.72)
    drawBoxSurfaceLocal(panel_limits, offset, [0.12 0.28 0.58], 0.52)
end

function drawBoxSurfaceLocal(limits, offset, face_color, face_alpha)
    x_min = limits(1,1) + offset(1);
    x_max = limits(1,2) + offset(1);
    y_min = limits(2,1) + offset(2);
    y_max = limits(2,2) + offset(2);
    z_min = limits(3,1) + offset(3);
    z_max = limits(3,2) + offset(3);

    vertices = [x_min y_min z_min;
                x_max y_min z_min;
                x_max y_max z_min;
                x_min y_max z_min;
                x_min y_min z_max;
                x_max y_min z_max;
                x_max y_max z_max;
                x_min y_max z_max];
    faces = [1 2 3 4;
             5 6 7 8;
             1 2 6 5;
             2 3 7 6;
             3 4 8 7;
             4 1 5 8];
    patch('Vertices',vertices,'Faces',faces,...
        'FaceColor',face_color,'FaceAlpha',face_alpha,...
        'EdgeColor',[0.08 0.10 0.12],'LineWidth',0.8)
end
