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

    sphLoad2 = load('GOESR_bus'); % Loading MSM model for GOES-R without boom

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
    
%% Plot all torques experienced by the servicer

PlotServicerTorqueCloud(params, relative_orientation_EMM_Torques);

%% Initial position plot

% This function does not update the params so using a different SN or DN is
% purely for making figures
EA = deg2rad([0,0,0]);
SN = Euler3212C(EA);
DN = eye(3);

[ ~, ~, ~, N_init_Lserv, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
    params.V, DN', SN' ,...
    params.debris.D_COM, params.servicer.S_COM);

N_L2_norm = N_init_Lserv/norm(N_init_Lserv);

PlotInitialPosition(params, SN, DN, [0,0,0]);

PlotInitialPositionContinuous(params, SN, DN, N_init_Lserv)

%% Anti_torque Attitude

[data_anti,i_anti,tot] = find_anti_torque(N_init_Lserv,relative_orientation_EMM_Torques);
        
[ ~, ~, ~, N_Lserv_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000, params.V,...
    DN', data_anti.SN', params.debris.D_COM, params.servicer.S_COM);

PlotAntiTorqueAttitude(params, DN, data_anti.SN, N_init_Lserv, N_Lserv_anti);

PlotAntiTorqueAttitudeContinuous(params, DN, data_anti.SN, N_init_Lserv, N_Lserv_anti)

PlotServicerInitialAntiOverlayContinuous(params, DN, SN, data_anti.SN, N_init_Lserv, N_Lserv_anti)

%% Reaction Wheel Simulations with E Torques

% Component of the Moment of Inertia of the reaction wheel about the wheel
% axis TODO: Justifiy wheel size
Iws = .12;

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
    tn = 1000*3600;
% Step size (s)
    dt = 1;
 
Ttot = 0:dt:tn;

params0 = params;

% Intial wheel speeds 
Om_0 = [0;0;0];

params.wheel_speed_threshold = 3000*(2*pi)/60;

% [data_anti,~,~] = find_anti_torque(H,data);

X0_servicer = [S_sig_BN0;S_w_BN;Om_0];
X0_target = [D_sig_BN0;D_w_BN];

% relative_orientation_EMM_Torques = [];
sim_flag = 0;
if sim_flag == 1
    params.wheel_speed_threshold = 6000*(2*pi)/60;
    
end



results =...
    N_RW_sim(X0_servicer,X0_target, Iws, Gs0, K, P, Ttot,params,relative_orientation_EMM_Torques);
%%
ShowPlots(params,results,Ttot(end),params.N_rvec_km,true)
