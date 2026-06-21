%% Electrostatic Momentem Management Base Code
clear;

% Add local code folders to the MATLAB path
addpath(genpath('MSM'))
addpath(genpath('functions'))

set(0,'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
% Load MSM models and Mass properties

% Parameters for SSL-1300 the acting target/debris
    params.debris.mass = 2000; % mass [kg]; 
    params.debris.voltage = -10e3; % [V] Voltage assuming fully conducting
    
    debris_D_COM = [0, 0, 0]'; % [m] SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
    
    debris_D_MI = [1000; 1000; 1000].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan

    sphLoad1 = load('SSL1300_bus.mat');% Loading MSM model for SSL-1300 geometry to match source link: body 2.8 x 2.1 x 2.0 m, panels 14 x 2.3 m each
    
% Parameters for GOES-R, the acting servicer/controlled spacecraft
    params.servicer.mass = 2857; % [kg] GOES-R series dry mass
    params.servicer.voltage = 10e3; % [V]
    
    % Parameters from the SolidWorks CAD in solid works frame
    SW_COM = [-2.174; 25.5; -31.469]*0.0254; % [m] COM 
    SW_MI = [39030563; 79464279; 94328208]*0.0002926396534292.*eye(3); % [kg m2] Moment of Inertia 

    % Maps from CAD axes to Solidworks axis ?
    DCM_IB = [0.014 0.965 -0.263; -0.041 0.263 0.964; 0.999 -0.003 0.043]; 

    % Maps from SoldiWorks model frame to SC body frame
    DCM_SW2B = [0 0 1; 1 0 0; 0 1 0]; 

    servicer_S_MI = DCM_SW2B*(DCM_IB*SW_MI*DCM_IB')*DCM_SW2B'; % in body frame
    servicer_S_COM = DCM_SW2B*SW_COM + [-2 0 0]'; % in body frame

    sphLoad2 = load('GOESR_bus_noboom'); % Loading MSM model for GOES-R without boom
%     sphLoad2 = load('GOESR_bus'); % Loading MSM model for GOES-R without boom


    % TODO: What are the offsets with the COM of mass?

% Proximity Parameters

% Intial Position with 0,0,0 at target COM
params.N_rvec_km = ([30 0 0]')./1000; % [km]

% DCM for a 90 deg rotation about the z axis 
M_90deg_Z = [cosd(90) , sind(90), 0;...
    -sind(90), cosd(90), 0;...
    0, 0, 1];

% DCM for a 90 deg rotation about the z axis 
M_90deg_X = [1,0,0;...
    0,cosd(90),-sind(90);...
    0,sind(90),cosd(90)];

% Rotation Matrix for a EA 2 rotation with EA = theta1
M_90deg_Y = [cosd(-90),0,-sind(-90);...
    0,1,0;...
    sind(-90),0,cosd(-90)];

DI = M_90deg_Y*M_90deg_X;

% Rotating the SSL craft to the intial orientation
for i = 1:length(sphLoad1.SPHSb)
    sph_loc = sphLoad1.SPHSb(1:3,i);
    new_sph_loc = DI*sph_loc;
    params.debris.D_spheres(1:3,i) = new_sph_loc;
    params.debris.D_spheres(4,i) = sphLoad1.SPHSb(4,i);
end
% params.debris.D_COM = debris_D_COM;
params.debris.D_COM = [2,0,0]';
params.debris.D_MI = debris_D_MI;

SI = M_90deg_Z*M_90deg_Y;

for i = 1:length(sphLoad2.SPHSb)
    sph_loc = sphLoad2.SPHSb(1:3,i);
    new_sph_loc = SI*sph_loc;
    params.servicer.S_spheres(1:3,i) = new_sph_loc;
    params.servicer.S_spheres(4,i) = sphLoad2.SPHSb(4,i);
end
% params.servicer.S_COM = servicer_S_COM;

params.servicer.S_COM = [2.7,-1,0]';

params.servicer.S_MI = servicer_S_MI;

% Vector of spacecraft potentials
params.V = [params.debris.voltage,params.servicer.voltage];

% Intial Rotation Matrix relative to the inital orientations.
params.servicer.NS = [1,0,0; 0,1,0;0,0,1];
params.debris.ND = [1,0,0; 0,1,0;0,0,1];


%% Additional Rotations

EA = deg2rad([-90,0,0]);
SN_i = Euler3212C(EA);
SN_i = eye(3);
DN_i = eye(3);

for i = 1:length(sphLoad2.SPHSb)
    sph_loc = sphLoad2.SPHSb(1:3,i);
    params.servicer.S_spheres(1:3,i) = SN_i'*params.servicer.S_spheres(1:3,i);
    params.servicer.S_spheres(4,i) = sphLoad2.SPHSb(4,i);
end

params.servicer.S_COM = SN_i'*params.servicer.S_COM;
params.servicer.S_MI = SN_i'*params.servicer.S_MI*SN_i'';


params.debris.D_COM = DN_i'*params.debris.D_COM;
params.debris.D_MI = DN_i'*params.debris.D_MI*DN_i'';


%% Computing Torques for all servicer Orientations (Stationary Target)

flag_solve4torques = 1;
n = 9;
clear relative_orientation_EMM_Torques
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000
if flag_solve4torques == 1
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    
elseif flag_solve4torques == 2
    params0 = params;
    % Debris modeled as a single 10 m radius sphere
    params.debris.mass = 2000; % [kg]
    params.debris.voltage = -10e3; % [V]
    params.debris.D_COM = [0;0;0]; % [m]
    params.debris.D_MI = 2/5*params.debris.mass*10^2*eye(3); % [kg m^2]
    params.debris.D_spheres = [0; 0; 0; 3];

    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
    params=params0;
elseif flag_solve4torques == 3
    % Debris modeled as five spheres aligned with the Y-axis
    params0 = params;
    params.debris.D_COM = [0;0;0]; % [m]
    Ix = 1/4*params.debris.mass*2.5^2 + 1/12*params.debris.mass*19.5^2;
    Iz = 1/4*params.debris.mass*2.5^2 + 1/12*params.debris.mass*19.5^2;
    Iy = .5*params.debris.mass*2.5^2;
    params.debris.D_MI = [Ix,0,0;0,Iy,0;0,0,Iz]; % [kg m^2]
    params.debris.D_spheres = [0 0 0 0 0;
        -14 -7 0 7 14;
        0 0 0 0 0;
        3 3 3 3 3];
    
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
    params=params0;
else
    load('n_50_30m_-10kV_SSL_10kV_GOESR_pos_1');
end
    
% Plot all torques experienced by the servicer

% PlotServicerTorqueCloud(params, relative_orientation_EMM_Torques);

%% Initial position plot


[ ~, ~, ~, N_init_Lserv, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
    params.V, eye(3), eye(3),...
    params.debris.D_COM, params.servicer.S_COM);

N_L2_norm = N_init_Lserv/norm(N_init_Lserv);

PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);

% PlotInitialPositionContinuous(params, SN, DN, N_init_Lserv)

% Computing the torques for all Target Rotations at first attitude

[rotating_debris_torques] = All_Torques_Target(params,n,eye(3));

PlotServicerTorqueCloud(params, rotating_debris_torques,SN_i);

% 
% torque_origin = params.servicer.S_COM + params.N_rvec_km*1000;
% figure
% hold on
% torque_norms = zeros(length(rotating_debris_torques),1);
% for i_L = 1:length(rotating_debris_torques)
%     torque_norms(i_L) = norm(rotating_debris_torques{i_L}.N_L_elect_serv);
% end
% max_torque = max(torque_norms);
% if max_torque == 0
%     max_torque = 1;
% end
% torque_arrow_scale = 8/max_torque;
% torque_colors = torqueMagnitudeColormapLocal(256);
% 
% for i_L = 1:length(rotating_debris_torques)
%     E_torque = rotating_debris_torques{i_L}.N_L_elect_serv;
%     if torque_norms(i_L) > 0
%         color_idx = min(256,max(1,ceil(256*torque_norms(i_L)/max_torque)));
%         E_torque = torque_arrow_scale*E_torque;
%         quiver3(torque_origin(1), torque_origin(2), torque_origin(3),...
%             E_torque(1),E_torque(2),E_torque(3),0,...
%             'Linewidth',0.8,'Color',torque_colors(color_idx,:))
%     end
% end
% 
% set(gca,'FontName','times')
% makeSphsPicture_2craft(params.debris.D_spheres, params.servicer.S_spheres,...
%     [0 0 0], params.N_rvec_km*1000,  [params.debris.voltage, params.servicer.voltage])
% 
% scatter3(params.debris.D_COM(1),params.debris.D_COM(2),params.debris.D_COM(3),20,'r','filled')
% scatter3(params.servicer.S_COM(1)+30,params.servicer.S_COM(2),params.servicer.S_COM(3),20,'k','filled')
% 
% % params.debris.D_COM, params.servicer.S_COM
% % c=colorbar;
% % c.Label.String = 'Surface Charge Density (nC/m^2)';
% % c.Label.FontSize = 14;
% 
% axis equal
% xlim([-3,40])
% ylim([-17,17])
% zlim([-12,12])
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% grid off
% view(3)
% hold off

%% Anti_torque Attitude

[data_anti,i_anti,tot] = find_anti_torque(N_init_Lserv,relative_orientation_EMM_Torques);
        
[ ~, ~, ~, N_Lserv_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000, params.V,...
    eye(3), data_anti.SN', params.debris.D_COM, params.servicer.S_COM);

PlotAntiTorqueAttitude(params, eye(3), data_anti.SN, N_init_Lserv, N_Lserv_anti);

% PlotAntiTorqueAttitudeContinuous(params, eye(3), data_anti.SN, N_init_Lserv, N_Lserv_anti)

% PlotServicerInitialAntiOverlayContinuous(params, eye(3), eye(3), data_anti.SN, N_init_Lserv, N_Lserv_anti)

%% Computing the torques for all Target Rotations

[rotating_debris_torques] = All_Torques_Target(params,n,data_anti.SN');

% PlotServicerTorqueCloud(params, rotating_debris_torques,data_anti.SN');

%% Reaction Wheel Simulations with E Torques

% Component of the Moment of Inertia of the reaction wheel about the wheel
% axis TODO: Justifiy wheel size
Iws = .2;

% Defining the "Gimbal" Frame
% Gimbal along the principle axes of the spacecraft
gs10 = [1;0;0];
gs20 = [0;1;0];
gs30 = [0;0;1];
Gs0 = [gs10,gs20,gs30];

% Issue with spacecraft voltage where you get the opposite of the expected
% results, i.e. opposite signs has repelling spacecraft

% Initial Servicer Orientation
S_sig_BN0 = [0,0,0]';
S_w_BN = [0;0;0];

% Initial Target Orientation
D_sig_BN0 = [0,0,0]';
D_w_BN = [0.0;0.0;0.0];


% Gains
    K = 5;
    P = 500;
% Total simulation time (s)
    tn = 10*3600;
% Step size (s)
    dt = 1;
 
Ttot = 0:dt:tn;

params0 = params;

% Intial wheel speeds 
Om_0 = [0;0;0];
Om_0 = [0;150;0];
params.sim.H = Om_0*Iws;
params.wheel_speed_threshold = 1000*(2*pi)/60;
% [data_anti,~,~] = find_anti_torque(H,data);

X0_servicer = [S_sig_BN0;S_w_BN;Om_0];
X0_target = [D_sig_BN0;D_w_BN];

% relative_orientation_EMM_Torques = [];
sim_flag = 0;
if sim_flag == 1
    params.wheel_speed_threshold = 6000*(2*pi)/60;
    
end

K_h = 0.0;

results =...
    N_RW_sim(X0_servicer,X0_target, Iws, Gs0, K, P, Ttot,params,relative_orientation_EMM_Torques,K_h);
%%
ShowPlots(params,results,Ttot(end),params.N_rvec_km,true)

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
