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
    params.debris.voltage = 10e3; % [V] Voltage assuming fully conducting
    
    debris_D_COM = [0, 0, 0]'; % [m] SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
    
    debris_D_MI = [1000; 1000; 1000].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan

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

    servicer_S_MI = DCM_SW2B*(DCM_IB*SW_MI*DCM_IB')*DCM_SW2B'; % in body frame
    servicer_S_COM = DCM_SW2B*SW_COM + [-2 0 0]'; % in body frame

    sphLoad2 = load('GOESR_bus_noboom'); % Loading MSM model for GOES-R without boom

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

% Computing Torques for all servicer Orientations (Stationary Target)

flag_solve4torques = true;

% 
% plotting_servicer = params.servicer.S_spheres;
%     for i = 1:size(params.servicer.S_spheres,2)
%         sph_loc = params.servicer.S_spheres(1:3,i);
%         N_sph_loc_serv = SN'*sph_loc;
%         plotting_servicer(1:3,i) = N_sph_loc_serv;
%     end
    
    
% Computing the Electrostatic Torques for different attitudes for SC 2
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000
if flag_solve4torques == true
    n = 9;
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
%     save('n_50_30m_-10kV_SSL_10kV_GOESR_pos_1','relative_orientation_EMM_Torques')
else
    load('n_50_30m_-10kV_SSL_10kV_GOESR_pos_1');
end

%% Initial position plot

EA = deg2rad([0,0,0]);
SN = Euler3212C(EA);
DN = eye(3);

% params.debris.N_spheres = [0,0,0,3]';
% params.servicer.N_spheres = [17.47-30,27.5-30;...
%     1,-5;...
%     1.5,-1.3;...
%     1,1];
% params.servicer.S_COM = [(17.47-30+27.5-30)/2,-2 ,0]';
% 
% params.debris.D_COM = [0,0,0]';
% 
% [N_F_on_debris, N_F_on_serv, N_L_elect_debris, N_L_elect_serv, qs, overlapFlag] = ...
%     multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
%     params.V, eye(3), eye(3), params.debris.N_COM, params.servicer.N_COM);

plotting_servicer = params.servicer.S_spheres;
        for i = 1:size(params.servicer.S_spheres,2)
            
            N_sph_loc_serv = SN'*params.servicer.S_spheres(1:3,i);
            plotting_servicer(1:3,i) = N_sph_loc_serv;
        end
               
        
[ ~, ~, ~, N_L2, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
    params.V, DN', SN' ,...
    params.debris.D_COM, params.servicer.S_COM);

serv_N_COM = SN'*params.servicer.S_COM;
deb_N_COM = DN'*params.debris.D_COM;


N_L2_norm = N_L2/norm(N_L2);

clf(figure(1))
figure(1)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.D_spheres, plotting_servicer,... %  plotting_servicer
[0 0 0], params.N_rvec_km*1000,  [params.debris.voltage, params.servicer.voltage])

quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, serv_N_COM(2) + params.N_rvec_km(2)*1000,serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_norm(1),N_L2_norm(2),N_L2_norm(3),6,'Linewidth',3)

scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')

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

% dots = [];
% for i = 1:length(relative_orientation_EMM_Torques)
% dots(i,1) = dot(N_L2/norm(N_L2),relative_orientation_EMM_Torques{i}.N_L_elect_serv/norm(relative_orientation_EMM_Torques{i}.N_L_elect_serv));
% end

[data_anti,i_anti,tot] = find_anti_torque(N_L2,relative_orientation_EMM_Torques);

% data_anti.SN = MRP2C([1,0,0]')

% Second position plot
plotting_servicer = params.servicer.S_spheres;
        for i = 1:size(params.servicer.S_spheres,2)
            plotting_servicer(1:3,i) = data_anti.SN'*params.servicer.S_spheres(1:3,i);
        end
        
        
[ ~, ~, ~, N_L2_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000, params.V,...
    DN', data_anti.SN', params.debris.D_COM, params.servicer.S_COM);
        
% data_anti.N_L_elect_serv = data_anti.SN*data_anti.N_L_elect_serv;
serv_N_COM = data_anti.SN'*params.servicer.S_COM;

N_L2_anti_norm = N_L2_anti/norm(N_L2_anti);

clf(figure(2))
figure(2)
hold on
set(gca,'FontName','times')
makeSphsPicture_2craft(params.debris.D_spheres, plotting_servicer,... % params.servicer.N_spheres plotting_servicer
[0 0 0], params.N_rvec_km*1000, [params.debris.voltage, params.servicer.voltage])

p(1)=quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, serv_N_COM(2) + params.N_rvec_km(2)*1000,serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_norm(1),N_L2_norm(2),N_L2_norm(3),6,'Linewidth',3);
p(2)=quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000,serv_N_COM(2) + params.N_rvec_km(2)*1000,serv_N_COM(3) + params.N_rvec_km(3)*1000,...
    N_L2_anti_norm(1),N_L2_anti_norm(2),N_L2_anti_norm(3),6,'Linewidth',3);

scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')

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
    tn = 100*3600;
% Step size (s)
    dt = 1;
 
Ttot = 0:dt:tn;

params0 = params;

% Intial wheel speeds 
Om_0 = [0;0;0];

params.wheel_speed_threshold = 500/60*(2*pi);

% [data_anti,~,~] = find_anti_torque(H,data);

X0_servicer = [S_sig_BN0;S_w_BN;Om_0];
X0_target = [D_sig_BN0;D_w_BN];

% relative_orientation_EMM_Torques = [];

results =...
    N_RW_sim(X0_servicer,X0_target, Iws, Gs0, K, P, Ttot,params,relative_orientation_EMM_Torques);
%%
ShowPlots(params,results,Ttot(end),params.N_rvec_km,true)
