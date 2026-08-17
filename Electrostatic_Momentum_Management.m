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

EA_s = deg2rad([0,0,0]);
SN_i = Euler3212C(EA_s);
% SN_i = eye(3);
EA_d = deg2rad([0,0,0]);
DN_i = Euler3212C(EA_d);
DN_i = eye(3);

for i = 1:length(sphLoad2.SPHSb)
    params.servicer.S_spheres(1:3,i) = SN_i'*params.servicer.S_spheres(1:3,i);
end

for i = 1:length(sphLoad1.SPHSb)
    params.debris.D_spheres(1:3,i) = DN_i'*params.debris.D_spheres(1:3,i);
end

params.servicer.S_COM = SN_i'*params.servicer.S_COM;
params.servicer.S_MI = SN_i'*params.servicer.S_MI*SN_i'';


params.debris.D_COM = DN_i'*params.debris.D_COM;
params.debris.D_MI = DN_i'*params.debris.D_MI*DN_i'';


%% Computing Torques for all servicer Orientations (Stationary Target)

flag_solve4torques = 1;
n = 9;
params0 = params;
clear relative_orientation_EMM_Torques
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000

torques_all_geo_T = [];
% for i_flag=1:4
% flag_solve4torques = i_flag;
if flag_solve4torques == 1
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
    save('All_torques_10kVs_30m_Sasym_Tsym.mat', 'relative_orientation_EMM_Torques')
elseif flag_solve4torques == 2

    % Debris modeled as a single 10 m radius sphere
    params.debris.mass = 2000; % [kg]
    params.debris.voltage = -10e3; % [V]
    params.debris.D_COM = [0;0;0]; % [m]
    params.debris.D_MI = 2/5*params.debris.mass*10^2*eye(3); % [kg m^2]
    params.debris.D_spheres = [0; 0; 0; 3];

    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
elseif flag_solve4torques == 3
    % Debris modeled as five spheres aligned with the Y-axis
    params.debris.D_COM = [0;0;0]; % [m]
    Ix = 1/4*params.debris.mass*3^2 + 1/12*params.debris.mass*34^2;
    Iz = 1/4*params.debris.mass*3^2 + 1/12*params.debris.mass*34^2;
    Iy = .5*params.debris.mass*3^2;
    params.debris.D_MI = [Ix,0,0;0,Iy,0;0,0,Iz]; % [kg m^2]
    params.debris.D_spheres = [1 1 1 1 1 1 1 1 1 1;
        -14.2000  -11.0444   -7.8889   -4.7333   -1.5778    1.5778    4.7333    7.8889   11.0444   14.2000;
        0 0 0 0 0 0 0 0 0 0;
        1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5];
    for i_sph = 1:10
        params.debris.D_spheres(1:3,i_sph) = DN_i'*params.debris.D_spheres(1:3,i_sph);
    end
    
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0],params0);
%     PlotInitialPositionContinuousCylinder(params, eye(3), eye(3), [0;0;0]);
elseif flag_solve4torques == 4
    params.debris.D_COM = DN_i'*[2.7,-1,0]';
    params.debris.D_MI = params.servicer.S_MI;
    params.debris.D_spheres = sphLoad2.SPHSb;
    for i = 1:length(sphLoad2.SPHSb)
        sph_loc = sphLoad2.SPHSb(1:3,i);
        new_sph_loc = DN_i'*SI*sph_loc;
        params.debris.D_spheres(1:3,i) = new_sph_loc;
        params.debris.D_spheres(4,i) = sphLoad2.SPHSb(4,i);
    end
    [relative_orientation_EMM_Torques, Ls,EAs] = All_Torques(params,n);
    PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
else
    load('All_torques_10kVs_30m_Sasym_Tsym.mat');
end
[ ~, ~, ~, Lserv, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
    params.V, eye(3), eye(3),...
    params.debris.D_COM, params.servicer.S_COM);

torques_all_geo_T = [torques_all_geo_T,Lserv];

  
%% Initial position plot

[Init_pos_targ_rot, ~, ~] = Avg_Servicer_Torque(params,n,eye(3));

N_L2_norm = Init_pos_targ_rot.Ls_target_rotation_norm;

PlotInitialPosition(params, eye(3), eye(3), N_L2_norm);

PlotServicerTorqueCloud(params, Init_pos_targ_rot.all_Ls_target_rotating, eye(3), flag_solve4torques);


%% Anti_torque Attitude

[data_anti,i_anti,tot] = find_anti_torque(N_L2_norm,relative_orientation_EMM_Torques);
        
[ ~, ~, ~, N_Lserv_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000, params.V,...
    eye(3), data_anti.SN', params.debris.D_COM, params.servicer.S_COM);

PlotAntiTorqueAttitude(params, eye(3), data_anti.SN, N_L2_norm, N_Lserv_anti);

[Anti_pos_targ_rot, ~, ~] = Avg_Servicer_Torque(params,n,data_anti.SN');


PlotServicerTorqueCloud(params, Anti_pos_targ_rot.all_Ls_target_rotating, data_anti.SN, flag_solve4torques);

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
    tn = 100*3600;
% Step size (s)
    dt = 1;
 
Ttot = 0:dt:tn;

params = params0;

% Intial wheel speeds 
Om_0 = [0;0;0];
% Om_0 = [0;150;0];
params.sim.H = Om_0*Iws;
params.wheel_speed_threshold = 3000*(2*pi)/60;
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




%%
% end
% %%
% serv_N_COM = params0.servicer.S_COM;
% deb_N_COM = params0.debris.D_COM;
% 
% figure
% hold on
% set(gca,'FontName','times')
% makeSphsPicture_2craft(params0.debris.D_spheres, params0.servicer.S_spheres,...
%     [0 0 0], params0.N_rvec_km*1000,  [params.debris.voltage, params.servicer.voltage])
% 
% 
% q(1) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,1),torques_all_geo_T(2,1),torques_all_geo_T(3,1),6e3,'Linewidth',3);
% 
% q(2) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,2),torques_all_geo_T(2,2),torques_all_geo_T(3,2),6e3,'Linewidth',3);
% 
% q(3) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,3),torques_all_geo_T(2,3),torques_all_geo_T(3,3),6e3,'Linewidth',3);
% 
% q(4) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,4),torques_all_geo_T(2,4),torques_all_geo_T(3,4),6e3,'Linewidth',3);
% 
% 
% scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
% scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')
% 
% 
% c=colorbar;
% c.Label.String = 'Surface Charge Density (nC/m^2)';
% c.Label.FontSize = 14;
% legend(q,'True Shape','Sphere', 'Cylinder','Single Panel')
% 
% axis equal
% xlim([27,40])
% ylim([-5,17])
% zlim([-12,12])
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% grid off
% view(3)
% hold off
% 
% figure
% hold on
% set(gca,'FontName','times')
% makeSphsPicture_2craft(params0.debris.D_spheres, params0.servicer.S_spheres,...
%     [0 0 0], params0.N_rvec_km*1000,  [params.debris.voltage, params.servicer.voltage])
% 
% 
% q(1) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,1),torques_all_geo_T(2,1),torques_all_geo_T(3,1),10e5,'Linewidth',3);
% 
% q(2) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,2),torques_all_geo_T(2,2),torques_all_geo_T(3,2),10e5,'Linewidth',3);
% 
% q(3) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,3),torques_all_geo_T(2,3),torques_all_geo_T(3,3),10e5,'Linewidth',3);
% 
% q(4) = quiver3(serv_N_COM(1) + params.N_rvec_km(1)*1000, ...
%     serv_N_COM(2) + params.N_rvec_km(2)*1000,...
%     serv_N_COM(3) + params.N_rvec_km(3)*1000,...
%     torques_all_geo_T(1,4),torques_all_geo_T(2,4),torques_all_geo_T(3,4),10e5,'Linewidth',3);
% 
% scatter3(deb_N_COM(1),deb_N_COM(2),deb_N_COM(3),20,'r','filled')
% scatter3(serv_N_COM(1)+30,serv_N_COM(2),serv_N_COM(3),20,'k','filled')
% 
% 
% c=colorbar;
% c.Label.String = 'Surface Charge Density (nC/m^2)';
% c.Label.FontSize = 14;
% legend(q,'Two Panel','Sphere', 'Cylinder','Single Panel')
% axis equal
% xlim([-10,40])
% ylim([-17,17])
% zlim([-12,12])
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
% grid off
% view(3)
% hold off
% 
%   



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
