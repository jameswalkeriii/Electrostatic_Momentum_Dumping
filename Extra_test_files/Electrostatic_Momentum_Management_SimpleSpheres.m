%% Electrostatic Momentum Management with Simple Sphere Geometry
clear;

% Debris modeled as a single 10 m radius sphere
params.debris.mass = 2000; % [kg]
params.debris.voltage = 10e3; % [V]
params.debris.D_COM = [0;0;0]; % [m]
params.debris.D_MI = 2/5*params.debris.mass*10^2*eye(3); % [kg m^2]
params.debris.N_spheres = [0; 0; 0; 3];

% Servicer modeled as three 3 m radius spheres aligned with the Y-axis
params.servicer.mass = 2857; % [kg]
params.servicer.voltage = -10e3; % [V]
params.servicer.S_COM = [0;6;0]; % [m]
Ix = 1/4*params.servicer.mass*2.5^2 + 1/12*params.servicer.mass*19.5^2;
Iz = 1/4*params.servicer.mass*2.5^2 + 1/12*params.servicer.mass*19.5^2;
Iy = .5*params.servicer.mass*2.5^2;
params.servicer.S_MI = [Ix,0,0;0,Iy,0;0,0,Iz]; % [kg m^2]
params.servicer.N_spheres = [0 0 0;
                             0 7 14;
                             0 0 0;
                             3 2.5 2.5];

% Relative position of the servicer origin in N
params.N_rvec_km = ([30 0 0]')./1000; % [km]

% Vector of spacecraft potentials
params.V = [params.debris.voltage, params.servicer.voltage];

% Initial DCMs relative to the initial orientations
params.servicer.SN = MRP2C([.5,-.2,.3]');
SN = params.servicer.SN;
params.debris.DN = eye(3);

%
flag_solve4torques = true;

% Computing the Electrostatic Torques for different attitudes for SC 2
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000
if flag_solve4torques == true
    n = 10;
    [relative_orientation_EMM_Torques, Ls,sigs,EAs] = All_Torques(params,n);
%     save('n_50_30m_-10kV_SSL_10kV_GOESR_pos_1','relative_orientation_EMM_Torques')
else
%     load('n_50_30m_-10kV_SSL_10kV_GOESR_pos_1');
end
% 
[ ~, ~, ~, N_L2, ~, overlapFlag] = ...
    multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.N_rvec_km*1000, params.V,...
    eye(3), SN, params.debris.D_COM, params.servicer.S_COM);
normL = N_L2/norm(N_L2);
% for i = 1:length(relative_orientation_EMM_Torques)
%     dots(i,1) = dot(normL,Ls(i,:));
% end

[data_anti,i_anti,tot] = find_anti_torque(N_L2,relative_orientation_EMM_Torques);

% Initial position plot
plotting_servicer = params.servicer.N_spheres;
        for i = 1:size(params.servicer.N_spheres,2)
            sph_loc = params.servicer.N_spheres(1:3,i);
            new_sph_loc_serv = SN*sph_loc;
            plotting_servicer(1:3,i) = new_sph_loc_serv;
        end

serv_N_COM = SN*params.servicer.S_COM;

% Initial position plot
clf(figure(1))
figure(1)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.N_spheres, plotting_servicer,...
[0 0 0], params.N_rvec_km*1000, ...
[params.debris.voltage, params.servicer.voltage])

scatter3(params.debris.D_COM(1),params.debris.D_COM(2),params.debris.D_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')

quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000,serv_N_COM(2) + params.N_rvec_km(2)*1000,serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2(1)*40000,N_L2(2)*40000,N_L2(3)*40000)
axis equal
xlim([-20,50])
ylim([-20,20])
zlim([-20,20])
c = colorbar;
c.Label.String = 'Surface Charge Density (nC/m^2)';
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
view(3)
hold off



% Second position plot
plotting_servicer = params.servicer.N_spheres;
        for i = 1:size(params.servicer.N_spheres,2)
            sph_loc = params.servicer.N_spheres(1:3,i);
            new_sph_loc_serv = data_anti.SN*sph_loc;
            plotting_servicer(1:3,i) = new_sph_loc_serv;
        end
        
        
[ ~, ~, ~, N_L2_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.N_rvec_km*1000, params.V,...
    eye(3), data_anti.SN,params.debris.D_COM, data_anti.SN*params.servicer.S_COM);
        
% data_anti.N_L_elect_serv = data_anti.C2*data_anti.N_L_elect_serv;
serv_N_COM = data_anti.SN*params.servicer.S_COM;

clf(figure(2))
figure(2)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.N_spheres, plotting_servicer,... % params.servicer.N_spheres plotting_servicer
[0 0 0], params.N_rvec_km*1000, [params.debris.voltage, params.servicer.voltage])

scatter3(params.debris.D_COM(1),params.debris.D_COM(2),params.debris.D_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')

quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000,serv_N_COM(2) + params.N_rvec_km(2)*1000,serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    data_anti.N_L_elect_serv(1)*40000,data_anti.N_L_elect_serv(2)*40000,data_anti.N_L_elect_serv(3)*40000)
axis equal
xlim([-3,50])
ylim([-17,17])
zlim([-12,12])
c=colorbar;
c.Label.String = 'Surface Charge Density (nC/m^2)';
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')

view(3)
hold off

%% Reaction Wheel Simulations with E Torques

Iws = .12;

gs10 = [1;0;0];
gs20 = [0;1;0];
gs30 = [0;0;1];
Gs0 = [gs10,gs20,gs30];

S_sig_BN0 = [0;0;0];
S_w_BN = [0;0;0];

D_sig_BN0 = [0;0;0];
D_w_BN = [0;0;0];

K = .0005;
P = .05;
tn = 4*3600;
dt = 1;
Ttot = 0:dt:tn;

Om_0 = [0;0;0];
params.wheel_speed_threshold = 500/60*(2*pi);

X0_servicer = [S_sig_BN0;S_w_BN;Om_0];
X0_target = [D_sig_BN0;D_w_BN];

relative_orientation_EMM_Torques = [];

results = ...
    N_RW_sim(X0_servicer, X0_target, Iws, Gs0, K, P, Ttot, params, relative_orientation_EMM_Torques);

ShowPlots(params, results, Ttot(end), params.N_rvec_km, true)
