function PlotServicerTorqueCloud(params, relative_orientation_EMM_Torques, C_plot, flag_solve4torques,SN_i)
% Plot all electrostatic torques experienced by the servicer.
% This function mirrors the inline plotting block in
% Electrostatic_Momentum_Management.m so both methods can be compared.

if nargin < 3
    C_plot = eye(3);
end
if nargin < 4
    flag_solve4torques = 1;
end

figure
hold on
set(gca,'FontName','times')

plotTargetGeometryLocal(params.debris.D_spheres, [0;0;0], eye(3), ...
    [0.72 0.76 0.80], 0.55, flag_solve4torques)
plotServicerSurfaceLocal(params.N_rvec_km*1000, C_plot, [0.16 0.32 0.66], 0.60)

torque_origin = C_plot'*SN_i*params.servicer.S_COM + params.N_rvec_km*1000;
scatter3(torque_origin(1),torque_origin(2),torque_origin(3),20,'k','filled')
debris_origin = params.debris.D_COM;
scatter3(debris_origin(1),debris_origin(2),debris_origin(3),20,'r','filled')

torque_norms = zeros(length(relative_orientation_EMM_Torques),1);
for i_L = 1:length(relative_orientation_EMM_Torques)
    torque_norms(i_L) = norm(relative_orientation_EMM_Torques{i_L}.N_L_elect_serv);
    Ls(i_L,:) = relative_orientation_EMM_Torques{i_L}.N_L_elect_serv;
end
max_torque = max(torque_norms);
if max_torque == 0
    max_torque = 1;
end
torque_arrow_scale = 10/max_torque;
torque_colors = torqueMagnitudeColormapLocal(256);

mean_torque_norm = torque_arrow_scale*mean(Ls);

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

quiver3(torque_origin(1), torque_origin(2), torque_origin(3),...
            mean_torque_norm(1),mean_torque_norm(2),mean_torque_norm(3),0,...
            'Linewidth',2,'Color',[0 0.4470 0.7410])

axis_limits = computePlotLimitsLocal(params, C_plot, torque_origin, ...
    relative_orientation_EMM_Torques, torque_arrow_scale, flag_solve4torques);
axis equal
xlim(axis_limits(1,:))
ylim(axis_limits(2,:))
zlim(axis_limits(3,:))
colormap(torque_colors)
caxis([0 max_torque])
cb = colorbar;
cb.Label.String = 'Torque Magnitude [N m]';
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid off
view(-45,20)
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

function plotTargetSurfaceLocal(SPHS, offset, DN, face_color, face_alpha)
    [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal(SPHS);
    drawBoxSurfaceLocal(bus_limits, offset, DN, face_color, face_alpha)
    drawBoxSurfaceLocal(panel_limits_1, offset, DN, face_color, face_alpha)
    drawBoxSurfaceLocal(panel_limits_2, offset, DN, face_color, face_alpha)
end

function plotTargetGeometryLocal(SPHS, offset, DN, face_color, face_alpha, flag_solve4torques)
    if flag_solve4torques == 3
        plotCylinderBodyLocal(SPHS, offset, DN, face_color, face_alpha)
    else
        plotTargetSurfaceLocal(SPHS, offset, DN, face_color, face_alpha)
    end
end

function plotCylinderBodyLocal(SPHS, offset, DN, face_color, face_alpha)
    [radius, y_limits, x_center, z_center] = cylinderParamsLocal(SPHS);

    n_theta = 40;
    n_y = 20;
    theta = linspace(0, 2*pi, n_theta);
    y_local = linspace(y_limits(1), y_limits(2), n_y);
    [Theta, Y] = meshgrid(theta, y_local);
    X = x_center + radius*cos(Theta);
    Yg = Y;
    Z = z_center + radius*sin(Theta);

    pts_B = [X(:)'; Yg(:)'; Z(:)'];
    pts_N = DN'*pts_B + offset(:);

    Xn = reshape(pts_N(1,:), size(X));
    Yn = reshape(pts_N(2,:), size(Yg));
    Zn = reshape(pts_N(3,:), size(Z));
    surf(Xn, Yn, Zn, 'FaceColor', face_color, 'FaceAlpha', face_alpha, ...
        'EdgeColor', 'none')

    cap_theta = linspace(0, 2*pi, n_theta);
    cap_r = linspace(0, radius, 16);
    [ThetaCap, RCap] = meshgrid(cap_theta, cap_r);
    for y_cap = y_limits
        Xc = x_center + RCap.*cos(ThetaCap);
        Yc = y_cap*ones(size(RCap));
        Zc = z_center + RCap.*sin(ThetaCap);
        cap_pts_B = [Xc(:)'; Yc(:)'; Zc(:)'];
        cap_pts_N = DN'*cap_pts_B + offset(:);
        Xcn = reshape(cap_pts_N(1,:), size(Xc));
        Ycn = reshape(cap_pts_N(2,:), size(Yc));
        Zcn = reshape(cap_pts_N(3,:), size(Zc));
        surf(Xcn, Ycn, Zcn, 'FaceColor', face_color, 'FaceAlpha', face_alpha, ...
            'EdgeColor', 'none')

        rim_pts_B = [x_center + radius*cos(cap_theta);
                     y_cap*ones(size(cap_theta));
                     z_center + radius*sin(cap_theta)];
        rim_pts_N = DN'*rim_pts_B + offset(:);
        plot3(rim_pts_N(1,:), rim_pts_N(2,:), rim_pts_N(3,:), ...
            'Color', [0.08 0.10 0.12], 'LineWidth', 0.8)
    end
end

function axis_limits = computePlotLimitsLocal(params, C_plot, torque_origin, ...
        relative_orientation_EMM_Torques, torque_arrow_scale, flag_solve4torques)
    % Collect the debris bus/panel bounds.
    debris_vertices = targetVerticesLocal(params.debris.D_spheres, eye(3), [0;0;0], ...
        flag_solve4torques);
    xyz_min = min(debris_vertices, [], 2);
    xyz_max = max(debris_vertices, [], 2);

    % Include debris COM marker.
    xyz_min = min([xyz_min, params.debris.D_COM], [], 2);
    xyz_max = max([xyz_max, params.debris.D_COM], [], 2);

    % Collect the rotated servicer bus/panel bounds.
    servicer_vertices = servicerVerticesLocal(C_plot, params.N_rvec_km*1000);
    xyz_min = min([xyz_min, min(servicer_vertices, [], 2)], [], 2);
    xyz_max = max([xyz_max, max(servicer_vertices, [], 2)], [], 2);

    % Include all torque arrow endpoints.
    xyz_min = min([xyz_min, torque_origin], [], 2);
    xyz_max = max([xyz_max, torque_origin], [], 2);
    for i_L = 1:length(relative_orientation_EMM_Torques)
        endpoint = torque_origin + torque_arrow_scale*relative_orientation_EMM_Torques{i_L}.N_L_elect_serv;
        xyz_min = min([xyz_min, endpoint], [], 2);
        xyz_max = max([xyz_max, endpoint], [], 2);
    end

    span = xyz_max - xyz_min;
    pad = max(0.1*span, 1);
    axis_limits = [xyz_min-pad, xyz_max+pad];
end

function vertices = servicerVerticesLocal(C_plot, offset)
    bus_limits = [0.55, 5.55;
                  -2.5, 2.5;
                  -2.5, 2.5];

    panel_limits = [1.05, 1.95;
                    -12.6, -2.5;
                    -2.4, 2.4];

    bus_vertices = boxVerticesLocal(bus_limits);
    panel_vertices = boxVerticesLocal(panel_limits);
    vertices = [bus_vertices, panel_vertices];
    vertices = C_plot'*vertices + offset(:);
end

function vertices = targetVerticesLocal(SPHS, DN, offset, flag_solve4torques)
    if flag_solve4torques == 3
        vertices = cylinderVerticesLocal(SPHS, DN, offset);
        return
    end

    [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal(SPHS);
    bus_vertices = boxVerticesLocal(bus_limits);
    panel_1_vertices = boxVerticesLocal(panel_limits_1);
    panel_2_vertices = boxVerticesLocal(panel_limits_2);
    vertices = [bus_vertices, panel_1_vertices, panel_2_vertices];
    vertices = DN'*vertices + offset(:);
end

function vertices = cylinderVerticesLocal(SPHS, DN, offset)
    [radius, y_limits, x_center, z_center] = cylinderParamsLocal(SPHS);
    vertices_B = [x_center - radius, y_limits(1), z_center - radius;
                  x_center + radius, y_limits(1), z_center - radius;
                  x_center + radius, y_limits(2), z_center - radius;
                  x_center - radius, y_limits(2), z_center - radius;
                  x_center - radius, y_limits(1), z_center + radius;
                  x_center + radius, y_limits(1), z_center + radius;
                  x_center + radius, y_limits(2), z_center + radius;
                  x_center - radius, y_limits(2), z_center + radius]';
    vertices = DN'*vertices_B + offset(:);
end

function [radius, y_limits, x_center, z_center] = cylinderParamsLocal(SPHS)
    x_center = mean(SPHS(1,:));
    z_center = mean(SPHS(3,:));
    radial_extent = sqrt((SPHS(1,:) - x_center).^2 + (SPHS(3,:) - z_center).^2) + SPHS(4,:);
    radius = max(radial_extent);
    y_limits = [min(SPHS(2,:) - SPHS(4,:)), max(SPHS(2,:) + SPHS(4,:))];
end

function vertices = boxVerticesLocal(limits)
    x_min = limits(1,1);
    x_max = limits(1,2);
    y_min = limits(2,1);
    y_max = limits(2,2);
    z_min = limits(3,1);
    z_max = limits(3,2);

    vertices = [x_min y_min z_min;
                x_max y_min z_min;
                x_max y_max z_min;
                x_min y_max z_min;
                x_min y_min z_max;
                x_max y_min z_max;
                x_max y_max z_max;
                x_min y_max z_max]';
end

function [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal(SPHS)
    %#ok<INUSD>
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
    vertices_B = boxVerticesLocal(limits);
    vertices_N = (BN'*vertices_B)' + offset(:)';
    faces = [1 2 3 4;
             5 6 7 8;
             1 2 6 5;
             2 3 7 6;
             3 4 8 7;
             4 1 5 8];
    patch('Vertices',vertices_N,'Faces',faces,...
        'FaceColor',face_color,'FaceAlpha',face_alpha,...
        'EdgeColor',[0.08 0.10 0.12],'LineWidth',0.8)
end
