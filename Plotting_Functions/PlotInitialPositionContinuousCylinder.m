function PlotInitialPositionContinuousCylinder(params, SN, DN, N_L2)
% Plot the servicer as a continuous bus/panel shape and the target as a
% continuous cylinder. The cylinder has length 31.4 m along the Y axis,
% radius 1.5 m, and local center at [1;0;0].

serv_N_COM = SN'*params.servicer.S_COM;
cyl_center_B = [1;0;0];
cyl_center_N = DN'*cyl_center_B;

if norm(N_L2) > 0
    N_L2_norm = N_L2/norm(N_L2);
else
    N_L2_norm = [0;0;0];
end

figure
hold on
set(gca,'FontName','times')

plotCylinderTargetLocal([0;0;0], DN, [0.72 0.76 0.80], 0.55)
plotCylinderOutlineLocal([0;0;0], DN, [0 0 0], 1.2)
plotTargetOutlineLocal([0;0;0], DN, [0 0 0], 1.4)
plotServicerSurfaceLocal(params.N_rvec_km*1000, SN, [0.16 0.32 0.66], 0.60)

quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
    serv_N_COM(2) + params.N_rvec_km(2)*1000, ...
    serv_N_COM(3) + params.N_rvec_km(3)*1000, ...
    N_L2_norm(1), N_L2_norm(2), N_L2_norm(3), 6, 'Linewidth', 3)

scatter3(cyl_center_N(1), cyl_center_N(2), cyl_center_N(3), 20, 'r', 'filled')
scatter3(serv_N_COM(1)+params.N_rvec_km(1)*1000, ...
    serv_N_COM(2)+params.N_rvec_km(2)*1000, ...
    serv_N_COM(3)+params.N_rvec_km(3)*1000, 20, 'k', 'filled')

axis equal
xlim([-5,40])
ylim([-20,20])
zlim([-7,7])
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
grid off
view(3)
camlight headlight
lighting gouraud
hold off
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

function plotCylinderTargetLocal(offset, DN, face_color, face_alpha)
radius = 1.5;
length_y = 31.4;
center_B = [1;0;0];

n_theta = 40;
n_y = 20;
theta = linspace(0, 2*pi, n_theta);
y_local = linspace(-length_y/2, length_y/2, n_y);
[Theta, Y] = meshgrid(theta, y_local);
X = center_B(1) + radius*cos(Theta);
Yg = center_B(2) + Y;
Z = center_B(3) + radius*sin(Theta);

pts_B = [X(:)'; Yg(:)'; Z(:)'];
pts_N = DN'*pts_B + offset(:);

Xn = reshape(pts_N(1,:), size(X));
Yn = reshape(pts_N(2,:), size(Yg));
Zn = reshape(pts_N(3,:), size(Z));
surf(Xn, Yn, Zn, 'FaceColor', face_color, 'FaceAlpha', face_alpha, ...
    'EdgeColor', 'none')

% End caps
cap_theta = linspace(0, 2*pi, n_theta);
cap_r = linspace(0, radius, 16);
[ThetaCap, RCap] = meshgrid(cap_theta, cap_r);
for y_cap = [-length_y/2, length_y/2]
    Xc = center_B(1) + RCap.*cos(ThetaCap);
    Yc = center_B(2) + y_cap*ones(size(RCap));
    Zc = center_B(3) + RCap.*sin(ThetaCap);
    cap_pts_B = [Xc(:)'; Yc(:)'; Zc(:)'];
    cap_pts_N = DN'*cap_pts_B + offset(:);
    Xcn = reshape(cap_pts_N(1,:), size(Xc));
    Ycn = reshape(cap_pts_N(2,:), size(Yc));
    Zcn = reshape(cap_pts_N(3,:), size(Zc));
    surf(Xcn, Ycn, Zcn, 'FaceColor', face_color, 'FaceAlpha', face_alpha, ...
        'EdgeColor', 'none')
end
end

function plotTargetOutlineLocal(offset, DN, edge_color, line_width)
    [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal();
    drawBoxOutlineLocal(bus_limits, offset, DN, edge_color, line_width)
    drawBoxOutlineLocal(panel_limits_1, offset, DN, edge_color, line_width)
    drawBoxOutlineLocal(panel_limits_2, offset, DN, edge_color, line_width)
end

function plotCylinderOutlineLocal(offset, DN, edge_color, line_width)
radius = 1.5;
length_y = 31.4;
center_B = [1;0;0];

theta = linspace(0, 2*pi, 160);
circle_x = center_B(1) + radius*cos(theta);
circle_z = center_B(3) + radius*sin(theta);

for y_cap = [-length_y/2, length_y/2]
    circle_pts_B = [circle_x;
                    center_B(2) + y_cap*ones(size(theta));
                    circle_z];
    circle_pts_N = DN'*circle_pts_B + offset(:);
    plot3(circle_pts_N(1,:), circle_pts_N(2,:), circle_pts_N(3,:), ...
        'Color', edge_color, 'LineWidth', line_width)
end
end

function [bus_limits, panel_limits_1, panel_limits_2] = targetSurfaceLimitsLocal()
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
patch('Vertices',vertices_N,'Faces',faces, ...
    'FaceColor',face_color,'FaceAlpha',face_alpha, ...
    'EdgeColor',[0.08 0.10 0.12],'LineWidth',0.8)
end

function drawBoxOutlineLocal(limits, offset, BN, edge_color, line_width)
vertices_B = boxVerticesLocal(limits);
vertices_N = (BN'*vertices_B')' + offset(:)';
edges = [1 2; 2 3; 3 4; 4 1;
         5 6; 6 7; 7 8; 8 5;
         1 5; 2 6; 3 7; 4 8];

for i_edge = 1:size(edges,1)
    pts = vertices_N(edges(i_edge,:), :);
    plot3(pts(:,1), pts(:,2), pts(:,3), '-', 'Color', edge_color, ...
        'LineWidth', line_width)
end
end

function vertices_B = boxVerticesLocal(limits)
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
end
