%% Electrostatic Momentem Management Base Code
clear;

% Set plot parameters
set(0,'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.5,'DefaultAxesTitleFontWeight', 'bold') ;

col_b = [0 0.4470 0.7410];
col_r = [0.8500 0.3250 0.0980];
col_y = [0.9290 0.6940 0.1250];
col_p = [0.4940 0.1840 0.5560];
col_g = [0.4660 0.6740 0.1880];
col_c = [0.3010 0.7450 0.9330];
col_m = [0.6350 0.0780 0.1840];
col_1 = [0.3350 0.5780 0.6840];
col_2 = [0.9350 0.4080 0.6840];
col_3 = [0.1350 0.2080 0.8840];
 

%% Load MSM models and Mass properties

% Parameters for SSL-1300 the acting target/debris
    params.debris.mass = 2000; % mass [kg]; 
    params.debris.voltage = 10e3; % [V] Voltage assuming fully conducting
    
    params.debris.COM = [0, 0, 0]'; % [m] SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
    
    params.debris.MI = [1000; 1000; 1000].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan

    sphLoad1 = load('SSL1300_bus.mat');% Loading MSM model for SSL-1300 geometry to match source link: body 2.8 x 2.1 x 2.0 m, panels 14 x 2.3 m each

% Parameters for GOES-R, the acting servicer/controlled spacecraft
    params.servicer.mass = 2857; % [kg] GOES-R series dry mass
    params.servicer.voltage = -10e3; % [V]
    
    % Parameters from the SolidWorks CAD in solid works frame
    SW_COM = [-2.174; 25.5; -31.469]*0.0254; % [m] COM 
    SW_MI = [39030563; 79464279; 94328208]*0.0002926396534292.*eye(3); % [kg m2] Moment of Inertia 

    % Maps from CAD axes to Solidworks axis ?
    DCM_IB = [0.014 0.965 -0.263; -0.041 0.263 0.964; 0.999 -0.003 0.043]; 

    % Maps from SoldiWorks model frame to SC body frame
    DCM_SW2B = [0 0 1; 1 0 0; 0 1 0]; 

    params.servicer.B_MI = DCM_SW2B*DCM_IB*SW_MI; % in body frame
    params.servicer.COM = DCM_SW2B*SW_COM + [-2 0 0]'; % in body frame

    sphLoad2 = load('GOESR_bus.mat'); % Loading MSM model for GOES-R

    % TODO: What are the offsets with the COM of mass?

% Proximity Parameters

% Setting Intial Position in the Hill Frame
initial_Xpos = 30; initial_Ypos = 0; initial_Zpos = 0; %[m]

% Intial Position with 0,0,0 at target COM
params.r_km = ([initial_Xpos initial_Ypos initial_Zpos]')./1000; % [km]

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

% Rotating the SSL craft to the intial orientation
for i = 1:length(sphLoad1.SPHSb)
    sph_loc = sphLoad1.SPHSb(1:3,i);
    new_sph_loc = M_90deg_Y*M_90deg_X*sph_loc;
    params.debris.spheres(1:3,i) = new_sph_loc;
    params.debris.spheres(4,i) = sphLoad1.SPHSb(4,i);
end

% Rotating the GOES-R craft to the intial orientation (3-1-2 EA rotation)

for i = 1:length(sphLoad2.SPHSb)
    sph_loc = sphLoad2.SPHSb(1:3,i);
    new_sph_loc = M_90deg_Z*M_90deg_Y*sph_loc;
    params.servicer.spheres(1:3,i) = new_sph_loc;
    params.servicer.spheres(4,i) = sphLoad2.SPHSb(4,i);
end

% Vector of spacecraft potentials
params.V = [params.debris.voltage,params.servicer.voltage];

% Intial Rotation Matrix relative to the inital orientations.
params.servicer.rotation_matrix = [1,0,0; 0,1,0;0,0,1];
params.debris.rotation_matrix = [1,0,0; 0,1,0;0,0,1];

%Initial position plot
clf(figure(1))
figure(1)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.servicer.spheres, params.debris.spheres, params.r_km*1000,...
[0 0 0], [params.servicer.voltage, params.debris.voltage])
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

%% Computing Torques for all servicer Orientations (Stationary Target)

flag_solve4torques = false;

% Computing the Electrostatic Torques for different attitudes for SC 2
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000
if flag_solve4torques == true
    n = 100;
    relative_orientation_EMM_Torques = All_Torques(params,n);
    save('n_100_30m_-10kV_SSL_10kV_GOESR_pos_1s','relative_orientation_EMM_Torques')
else
    load('n_100_30m_-10kV_SSL_10kV_GOESR_pos_1');
end

[ ~, ~, ~, L2, ~, overlapFlag] = ...
    multisphereFT( params.debris.spheres, params.servicer.spheres, params.r_km*1000, params.V,...
    eye(3), eye(3),params.debris.COM, params.servicer.COM);

[data_anti,i_anti,tot] = find_anti_torque(L2,relative_orientation_EMM_Torques);
 
%% Reaction Wheel Simulations with E Torques

% Component of the Moment of Inertia of the reaction wheel about the wheel
% axis TODO: Justifiy wheel size
Iws = .12;

% The moment of Inertia of the spacecraft excluding the reaction wheels
I_RW = params.servicer.B_MI;

% Defining the "Gimbal" Frame
% Gimbal along the principle axes of the spacecraft
gs10 = [1;0;0];
gs20 = [0;1;0];
gs30 = [0;0;1];
Gs0 = [gs10,gs20,gs30];

% Issue with spacecraft voltage where you get the opposite of the expected
% results, i.e. opposite signs has repelling spacecraft

% Distance between center of masses is 8 m: spacecraft can collide 

% Number of Reaction Wheels
N = 3;

% Intial Servicer Orientation
B_sig_BN0 = [0,0,0]';
B_w_BN = [0;0;0];


% Gains
    K = 5;
    P = 500;
% Total simulation time (s)
    tn = 12*3600;
% Step size (s)
    dt = 1;

params0 = params;


% Intial wheel speeds 
Om_0 = [0;0;0];

params.wheel_speed_threshold = 5000/60/(2*pi);

% [data_anti,~,~] = find_anti_torque(H,data);

X0 = [B_sig_BN0;B_w_BN;Om_0];


[Xtot, Htot, Ttot, ttot, aterrtot, werrtot, ustot, Lrtot,L_e_tot,params,...
    normBN,reftot,rtot,flags] =...
    N_RW_sim(X0, I_RW, Iws, Gs0, N, K, P, tn, dt,params,relative_orientation_EMM_Torques);
%%
ShowPlots(Xtot, Htot, Ttot, ttot, aterrtot, werrtot,...
    ustot,Lrtot,L_e_tot,params0,params,normBN,reftot,tn,params.r_km,true)


