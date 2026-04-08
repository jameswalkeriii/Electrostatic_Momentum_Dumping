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
params.servicer.SN = eye(3);
params.debris.DN = eye(3);

% Initial position plot
clf(figure(1))
figure(1)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.N_spheres, params.servicer.N_spheres,...
[0 0 0], params.N_rvec_km*1000, ...
[params.debris.voltage, params.servicer.voltage])
scatter3(params.debris.D_COM(1),params.debris.D_COM(2),params.debris.D_COM(3),20,'r','filled')
scatter3(params.servicer.S_COM(1)+30,params.servicer.S_COM(2),params.servicer.S_COM(3),20,'k','filled')
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
