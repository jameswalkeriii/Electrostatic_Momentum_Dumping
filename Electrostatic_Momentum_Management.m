%% Electrostatic Momentem Management Base Code
clear;
Colors_for_Plotting
% Add local code folders to the MATLAB path
addpath(genpath('MSM'))
addpath(genpath('attitude_functions'))
addpath(genpath('Data_Bases'))
addpath(genpath('Plotting_Functions'))
addpath(genpath('support'))

set(0,'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)

%% Additional Rotations - Changing First Attitude
EA_s = deg2rad([-90,0,0]); SN_i = Euler3212C(EA_s);
EA_d = deg2rad([0,0,0]); DN_i = Euler3212C(EA_d);

params.NS_i = SN_i';
params.ND_i = DN_i';

%% Loading and Rotating Servicer object to default position
servicer_shape = "GOESR_bus_noboom";
params.servicer = load_Servicer(servicer_shape,params.NS_i);

%% Loading and Rotating Debris object to default position
derbis_shape = "SSL1300_bus";
params.debris = load_Debris(derbis_shape,params.ND_i);

%% Charging and Spacecraft Scenario Parameters

% Intial seperation vector - defines from debris COM to servicer COM
params.N_rvec_km = ([30 0 0]')./1000; % [km]

% Spacecraft Potentials assuming fully conducting spacecraft
params.servicer.voltage = 10 *1e3; % [V]
params.debris.voltage = -10 *1e3; % [V] 

% Vector of spacecraft potentials
params.V = [params.debris.voltage,params.servicer.voltage];

PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);
%% Computing Average Torques for all servicer for a servicer at all attiudes

% n has to be an even number so -180 - +180 is n steps and -90 - +90 is n/2 steps potentially have more target attitudes than servier attitudes
    % Defined using 3-2-1 EAs
flag_build_database = 1; % Used for defining the torque database 
                            % 0 - Used stored data based
                            % 1 - Build torque data base

% TODO: Fix sphere and cylinder shapes so they
% Current Build shapes do not rotate with DN_i
build_shape = params.debris.shape; % Defining the target shape used to build torque database
                                        % params.debris.shape - perfect knowledge of the target
                                        % "Sphere" 
                                        % "Cylinder" 
 
% n attitudes - gives data base of n^3/2
n = 6;

% Store default parameters
params0 = params;

if flag_build_database == 0
    switch build_shape
        case "SSL1300_bus"
            % TODO: Genearilze name to spacecraft seperation distance,potential, and shapes
            load("All_torques_10kVs_30m_Sasym_Tsym_n" +n+".mat");
                  
        case "GOESR_bus_noboom"
            
        case "Sphere"
            load("All_torques_10kVs_30m_Sasym_3mSphere_n" +n+".mat");
        case "Cylinder"
            load("All_torques_10kVs_30m_Sasym_31p4mCyl_n" +n+".mat");
    end
elseif flag_build_database == 1
    
    params.debris = load_Debris(build_shape,params.ND_i);
    
    % Intial seperation vector - defines from debris COM to servicer COM
    params.N_rvec_km = ([30 0 0]')./1000; % [km]
    
    % Spacecraft Potentials assuming fully conducting spacecraft
    params.servicer.voltage = 10 *1e3; % [V]
    params.debris.voltage = -10 *1e3; % [V]
    
    % Vector of spacecraft potentials
    params.V = [params.debris.voltage,params.servicer.voltage];
    
    % Save data so only needs to be computed once
    switch build_shape
        case "SSL1300_bus"
            [relative_orientation_EMM_Torques, ~,~] = All_Torques(params,n);
%             save("All_torques_10kVs_30m_Sasym_Tsym_n"+n+".mat", 'relative_orientation_EMM_Torques', '-v7.3');
        case "GOESR_bus_noboom"
            [relative_orientation_EMM_Torques, ~,~] = All_Torques(params,n);
        case "Sphere"
            [relative_orientation_EMM_Torques, ~,~] = All_Torques(params,n);
            save("All_torques_10kVs_30m_Sasym_3mSphere_n" +n+".mat", 'relative_orientation_EMM_Torques', '-v7.3');
        case "Cylinder"   
            [relative_orientation_EMM_Torques, ~,~] = All_Torques(params,n);
            save("All_torques_10kVs_30m_Sasym_31p4mCyl_n" +n+".mat", 'relative_orientation_EMM_Torques', '-v7.3');
    end
end
       
    
% Plot Intial positions of PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);    
PlotInitialPosition(params, eye(3), eye(3), [0,0,0]);

% Reset debris parameters back to default now that the database has been built
params = params0;
ada = [];
sd = [];
for i=1:length(relative_orientation_EMM_Torques)
    ada(i) = relative_orientation_EMM_Torques{i}.average_dist_angle;
    sd(i) = relative_orientation_EMM_Torques{i}.std_dist_angle;
    
end
%% Removing "Bad" attitudes from the database

% iadd = 1;
% for i_all = 1:n^2*round(n/2)
%     if relative_orientation_EMM_Torques{i_all}.dist == 1
%         bad_atts{iadd} = relative_orientation_EMM_Torques{i_all};
%         relative_orientation_EMM_Torques{i_all} = {};
%         iadd = iadd + 1;
%     end
% end

relative_orientation_EMM_Torques = relative_orientation_EMM_Torques(~cellfun('isempty', relative_orientation_EMM_Torques));

%% Rotate Database torques to be in the now rotated servicer frame (if NS_i =! 0)
% R = params.NS_i;
% % R = eye(3);
% for i = 1:length(relative_orientation_EMM_Torques)
%     relative_orientation_EMM_Torques{i}.avg_N_L_elect_serv = R * relative_orientation_EMM_Torques{i}.avg_N_L_elect_serv(:);
%     
%     relative_orientation_EMM_Torques{i}.NS = relative_orientation_EMM_Torques{i}.NS*R';
%     relative_orientation_EMM_Torques{i}.SN = relative_orientation_EMM_Torques{i}.NS;
% end

%% Initial position plot

[Init_pos_targ_rot, ~, ~,Init_pos_targ_rot_all_Ls_target_rotating] = Avg_Servicer_Torque(params,n,eye(3));

N_L2_norm = Init_pos_targ_rot.Ls_target_rotation_norm;

PlotInitialPosition(params, eye(3), eye(3), N_L2_norm);

% Fix/ensure PlotServicerTorqueCloud works properly
% PlotServicerTorqueCloud(params, Init_pos_targ_rot_all_Ls_target_rotating, eye(3), flag_build_database);

%% Anti_torque Attitude

[data_anti,i_anti,tot] = find_anti_torque(N_L2_norm,relative_orientation_EMM_Torques,params.NS_i);
        
[ ~, ~, ~, N_Lserv_anti, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000, params.V,...
    eye(3), data_anti.SN', params.debris.D_COM, params.servicer.S_COM);

PlotAntiTorqueAttitude(params, eye(3), data_anti.SN, N_L2_norm, N_Lserv_anti);

[Anti_pos_targ_rot, ~, ~,Anti_pos_targ_rot_all_Ls_target_rotating] = Avg_Servicer_Torque(params,n,data_anti.SN');


% PlotServicerTorqueCloud(params, Anti_pos_targ_rot_all_Ls_target_rotating, data_anti.SN, flag_build_database);

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
    K = 2;
    P = 1500;
% Total simulation time (s)
    tn = 600*3600;
% Step size (s)
    dt = 2;
 
Ttot = 0:dt:tn;

% Intial wheel speeds 
Om_0 = [0;0;0];
% Om_0 = N_L2_norm'*650;

params.sim.H = Om_0*Iws;
params.wheel_speed_threshold = 3000*(2*pi)/60;

X0_servicer = [S_sig_BN0;S_w_BN;Om_0];
X0_target = [D_sig_BN0;D_w_BN];

% relative_orientation_EMM_Torques = [];
sim_flag = 0;

K_h = 0.0;

results =...
    N_RW_sim(X0_servicer,X0_target, Iws, Gs0, K, P, Ttot,params,relative_orientation_EMM_Torques,K_h);

ShowPlots(params,results,Ttot(end),params.N_rvec_km,true)
