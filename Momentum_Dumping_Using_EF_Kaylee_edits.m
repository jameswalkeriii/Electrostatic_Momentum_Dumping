 %% Momentum Dumping using Electrostatic Forces
clear;
set(0,'defaulttextinterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex')
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
set(0,'DefaultAxesTitleFontSizeMultiplier', 1.5,'DefaultAxesTitleFontWeight', 'bold') ;
 

%% Load MSM models and Mass properties
% Estimated COM positions for each, relative to MSM model origin

% SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
params.COM1 = [0, 0, 0]'; 
%GOES-R ish COM, in inches converted to meters
params.COM2 = [-2.174; 25.5; -31.469]*0.0254; 

% SSL-1300 mass; geometry to match source link: body 2.8 x 2.1 x 2.0 m, panels 14 x 2.3 m each
params.m1 = 2000; 
% kg, GOES-R series dry mass
params.m2 = 2857; 

% Moment of Inertia found from solidworks, in kg-m^2
params.I2 = [39030563;
    79464279;
    94328208]*0.0002926396534292.*eye(3);

% Maps from CAD axes to body axis
DCM_IB = [0.014 0.965 -0.263;
    -0.041 0.263 0.964;
    0.999 -0.003 0.043]; 

% Inertia in SC body axes
params.I2 = DCM_IB*params.I2;

% Maps from SoldiWorks model frame to SC body frame
DCM_SW2B = [0 0 1;
    1 0 0;
    0 1 0]; 

params.I2 = DCM_SW2B*params.I2; %in body frame
params.COM2 = DCM_SW2B*params.COM2 + [-2 0 0]';

% SSL Moment of Inertia From email with Dan
params.I1 = [1000; 1000; 1000].*eye(3); 

% Loading MSM model for SSlL-1300
sphLoad1 = load('SSL1300_bus.mat');
params.SPHS1 = sphLoad1.SPHSb;

% Loading MSM model for GOES-R
sphLoad2 = load('GOESR_bus.mat');
params.SPHS2 = sphLoad2.SPHSb;

% Setting Intial Position in the Hill Frame
initialXpos = 30; % [m]
initialYpos = 0;
initialZpos = 0;

% Hill Frame Intial Position
x0_H = ([initialXpos initialYpos initialZpos]')./1000;

%Craft 1 voltage
params.V1 = 10e3;
%Craft 2 voltage
params.V2 = -10e3; 

% Defining the Euler Angles that define the rotation of the MSM models into
% their desired intial orientation
theta1 = deg2rad(90);
theta2 = deg2rad(0);
theta3 = deg2rad(90);

% Rotation Matrix for a EA 3 rotation with EA = theta1
M3_theta1 = [cos(theta3) , sin(theta3), 0;...
    -sin(theta3), cos(theta3), 0;...
    0, 0, 1];

% Rotation Matrix for a EA 1 rotation with EA = theta1
M1_theta1 = [1,0,0;...
    0,cos(theta1),-sin(theta1);...
    0,sin(theta1),cos(theta1)];

% Rotation Matrix for a EA 2 rotation with EA = theta1
M2_theta1 = [cos(-theta2),0,-sin(-theta2);...
    0,1,0;...
    sin(-theta2),0,cos(-theta2)];

% Rotation Matrix for a EA 3 rotation with EA = theta4
M3_theta4 = [cos(pi/1.01) , sin(pi/1.01), 0;...
    -sin(pi/1.01), cos(pi/1.01), 0;...
    0, 0, 1];

% Rotation Matrix for a EA 1 rotation with EA = theta4
M1_theta4 = [1,0,0;...
    0,cos(pi/4),-sin(pi/4);...
    0,sin(pi/4),cos(pi/4)];

% DCM_negL =  M2_theta1*M1_theta1*M3_theta1;
M1_90 = [1,0,0;...
    0,cos(180),-sin(180);...
    0,sin(180),cos(180)];
% Rotating the SSL craft to the intial orientation
for i = 1:length(params.SPHS1)
    sph_loc = params.SPHS1(1:3,i);
    new_sph_loc = sph_loc;
    params.SPHS1(1:3,i) = new_sph_loc;
end

% Rotating the GOES-R craft to the intial orientation (3-1-2 EA rotation)
% C3 = MRP2C([1/3,-1/3,1/3]);
for i = 1:length(params.SPHS2)
    sph_loc = params.SPHS2(1:3,i);
    new_sph_loc =MRP2C([-1;0;0])*sph_loc;
    params.SPHS2(1:3,i) = new_sph_loc;
end

% Intial Separation Vector
params.r = [initialXpos initialYpos initialZpos]';
% Vector of spacecraft potentials
params.V = [params.V1,params.V2];
% Intial Rotation Matrix relative to the inital orientations.
params.C1 = [1,0,0; 0,1,0;0,0,1];

params.C2 = [1,0,0; 0,1,0;0,0,1];

% % Calculate Initial Electrostatic Forces (F1,  F2) and Torqus (L1, L2) on
% % Spacecraft 1 and 2 respectively. qs is the vector of charge on all 
% % spheres. OverlapFlag is if any of the spheres overlap
% [ F1, F2, L1, L2, qs, overlapFlag] = ...
%     multisphereFT( params.SPHS1, params.SPHS2, params.r, params.V,...
%     params.C1, params.C2,params.COM1, params.COM2);
%
% %Initial position plot
% figure 
% grid on
% hold on
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, params.SPHS2, x0_H*1000,...
% [0 0 0], params.V1)
% hold on
% axis tight
% c=colorbar;
% c.Label.String = 'Surface Charge Density (nC/m^2)';
% quiver3(0,0,0,15^3*L2(1),15^3*L2(2),15^3*L2(3),'k','Linewidth',2)
% axis tight
% % title('Initial Positions')
% 
% view(3)
% hold off
%
% % Rotating the Inital Position of SC2 by 180 degress
% M3_180 = [cosd(180) , sind(180), 0;...
%     -sind(180), cosd(180), 0;...
%     0, 0, 1];
% 
% M2_180 = [cosd(180),0,-sind(180);...
%     0,1,0;...
%     sind(180),0,cosd(180)];
% M1_180 = [1,0,0;...
%     0,cosd(180),-sind(180);...
%     0,sind(180),cosd(180)];
% 
% inter_sph_loc = params.SPHS2*0;
% for i = 1:length(params.SPHS2)
%     inter_sph_loc(1:3,i) = M2_180*M3_180*params.SPHS2(1:3,i);
%     inter_sph_loc(4,i) = params.SPHS2(4,i);
% end

% [ ~, ~, L1_M3, L2_M3, qs, overlapFlag] = ...
%     multisphereFT( params.SPHS1, inter_sph_loc, params.r, params.V,...
%     params.C1, params.C2,params.COM1, params.COM1);

% %Initial position plot with 180 rotation
% figure 
% grid on
% hold on
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, inter_sph_loc, x0_H*1000,...
% [0 0 0], params.V1)
% hold on
% quiver3(0,0,0,10000*L2_M3(1),10000*L2_M3(2),10000*L2_M3(3),'k','Linewidth',2)
% title('Initial Positions with 180 rotation')
% view(3)
% axis tight
% 
% hold off

% look into papers that talk about movie spinning 3D shapes is easier than
% a slowly rotating one 


%%
% Finding the attitude with the Electrostatic Torque most opposite to the
% Electrostatic Torque, L2, computed at the intial, unrotated position.

% Does not currently work
% L2 = -[0.0003, -0.0026, 0.0018]';
% [data_anti,i_anti,tot] = find_anti_torque(L2,data);
% %%
% Rotating the spheres of the second spacecraft to represent the
% anti-torque attitude
% C2_anti = data_anti.C2;
% anti_sph_loc = params.SPHS2*0;
% for i = 1:length(params.SPHS2)
%     anti_sph_loc(1:3,i) = C2_anti*params.SPHS2(1:3,i);
%     anti_sph_loc(4,i) = params.SPHS2(4,i);
% end
% L_anti = data_anti.T_L;
% % Computing the Electrostatic Torques at the computed anti-torque attitude
% % At the current anti-torque positon ([-0.8422  0.5046  -0.1001])
% [ ~, ~, L1_a, L2_a, qs, overlapFlag_a] = ...
%     multisphereFT( params.SPHS1, params.SPHS2, params.r, params.V,...
%     params.C1, C2_anti,params.COM1, params.COM1);
% 
% % Servicing attitude
% figure 
% grid on
% hold on
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, params.SPHS2, x0_H*1000,...
% [0 0 0], params.V1)
% hold on
% h = quiver3(0,0,0,15^3*L2(1),15^3*L2(2),15^3*L2(3),'r','Linewidth',2,'MaxHeadSize',.8);
% t = hgtransform('Parent',gca);
% R = makehgtform('axisrotate',[L2(1),L2(2),L2(3)],pi/4);
% set(t,'Matrix',R);
% set(h,'Parent',t); % The arrows should point correctly now
% % title('Initial Positions with 180 rotation')
% view(3)
% axis tight
% hold off
% 
% % Plotting the "anti-torque" attitude using MSM
% figure 
% grid on
% hold on
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, anti_sph_loc, x0_H*1000,...
% [0 0 0], params.V1)
% hold on
% h = quiver3(0,0,0,15^3*L_anti(1),15^3*L_anti(2),15^3*L_anti(3),'b','Linewidth',2,'MaxHeadSize',.8);
% t = hgtransform('Parent',gca);
% R = makehgtform('axisrotate',[L2_M3(1),L2_M3(2),L2_M3(3)],pi/4);
% set(t,'Matrix',R);
% set(h,'Parent',t); % The arrows should point correctly now
% 
% % title('Anti L Position')
% % c=colorbar;
% % c.Label.String = 'Surface Charge Density (nC/m^2)';
% % c.Label.FontSize = 14;
% view(3)
% axis tight
% hold off

 
%% Reaction Wheel Simulations with E Torques

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

% (1,1) Component of the Moment of Inertia of the reaction wheel about the
% wheel axis 
Iws = .12;
% Iws = 1;
% The moment of Inertia of the spacecraft excluding the reaction wheels
I_RW = params.I2;

% Defining the "Gimbal" Frame
% Gimbal along the principle axes of the spacecraft
gs10 = [1;0;0];
gs20 = [0;1;0];
gs30 = [0;0;1];
Gs0 = [gs10,gs20,gs30];



% Spacecraft 1 Voltage
params.V1 = 10e3;
% Spacecraft 2 Voltage
params.V2 = -10e3;
% Issue with spacecraft voltage where you get the opposite of the expected
% results, i.e. opposite signs has repelling spacecraft
params.V = [params.V1,params.V2];

% Distance between center of masses is 8 m: spacecraft can collide 
params.r = [30;0;0];

% Number of Reaction Wheels
N = 3;

% Selelct Intial Servicer Orientations

% theta4 = -pi/3;
% M3_theta4 = [cos(theta4) , sin(theta4), 0;...
%     -sin(theta4), cos(theta4), 0;...
%     0, 0, 1];
% B_sig_BN0 = C2MRP(M3_theta4*M2_theta1*M1_theta1*M3_theta1)';
% B_sig_BN0 = [-.5142, -.213,.5142]';
% B_sig_BN0 = [1/3,-1/3,1/3]';
% B_sig_BN0 = [0.5142,-.2130,-.5142]';
% B_sig_BN0 = [.3403,-.3298,-3403]';
B_sig_BN0 = [0,0,0]';
% B_sig_BN0 = [1,0,0]';

% Initial Angular Velocity of the servicer
B_w_BN = [0;0;0];

% params.SPHS2 = inter_sph_loc;

% State Vector


% Gains
K = 5;
P = 500;
% Total simulation time (s)
% ttotal = 200000;
ttotal = 10*24*3600;
tn = ttotal;
% Step size (s)
dt = 10;

params0 = params;
% Om_0 = -[280.4169;208.9686;3.5297];
% % Initial Angular Momentum
% HB = [0;0;0];
% 
% HW = Iws*(Om_0(1))*gs10 + Iws*(Om_0(2))*gs20 + Iws*(Om_0(3))*gs30;
% H = HB + HW;

% Intial wheel speeds 
% Wheels Intially not at rest
% Om_0 = [89.4726;168.4593;-17.4980];
% Om_0 = -[280.4169/1.6;208.9686/1.1;3.5297/.9];
% Om_0 = -[280.4169;208.9686;3.5297];
Om_0 = [0;0;0];

% [data_anti,~,~] = find_anti_torque(H,data);
% 
% % sig_RN = C2MRP(data_anti.C2);
% % sig_RN = [-0.4731, -0.5171, -0.4917]'; %8777
% sig_RN = C2MRP(data{68567}.C2); %68567
% B_sig_BN0 = sig_RN;
% Om_0 = -[280.4169;208.9686;3.5297];
X0 = [B_sig_BN0;B_w_BN;Om_0];
% [Xtot, Htot, Ttot, ttot,aterrtot,werrtot,ustot,Lrtot,...
%     L_e_tot,params,normBN,reftot,rtot,flags] = N_RW_sim_no_man(X0, I_RW, Iws, Gs0,...
%     N, K, P, tn, dt,params,data,sig_RN);
% % Show Plots
% flag_plots = true;
% % Make Video: takes a really long time, only use when necessary
% flag_vid = false;
% % Remember to change file name for each video
% filenam = 'Etorques_Tumble_testing_opposite_Torque';
% %
% if flag_plots == true
%     ShowPlots(Xtot, Htot, Ttot, ttot, aterrtot, werrtot,ustot,Lrtot,L_e_tot,...
%     params0,params,normBN,reftot,ttotal,x0_H,true) 
% end
% if flag_vid == true
%     MakeVideo(Xtot,params,x0_H,filenam,rtot)
% end
% 
% time_mat = [10;...
%     54140];

[ ~, ~, L1_1, L2_1, qs, overlapFlag] = ...
    multisphereFT( params.SPHS1, params.SPHS2, params.r, params.V,...
    params.C1, params.C2,params.COM1, params.COM2);
nscale = 10000;
% Plotting the Initial attitude using MSM
figure 
grid on
hold on
xlabel('Relative X distance [m]')
ylabel('Relative Y distance [m]')
zlabel('Relative Z distance [m]')
set(gca,'FontName','times')
makeSphsPicture_2craft( params.SPHS1, params.SPHS2, x0_H*1000,...
[0 0 0], [params.V1,params.V1]')
hold on
quiver3(0,0,0,L2_1(1,1)*nscale,L2_1(2,1)*nscale,L2_1(3,1)*nscale,'color',col_b,'Linewidth',3,'MaxHeadSize',1)
% c=colorbar;
% c.Label.String = 'Surface Charge Density (nC/m^2)';
% c.Label.FontSize = 14;
view(3)
axis tight
hold off



SPHS2_2(1:3,:) = MRP2C([-1;0;0])*params.SPHS2(1:3,:);
SPHS2_2(4,:) = params.SPHS2(4,:);

[ ~, ~, L1_2, L2_2, qs, overlapFlag] = ...
    multisphereFT( params.SPHS1, SPHS2_2, params.r, params.V,...
    params.C1, params.C2,params.COM1, params.COM2);

% Plotting the Second attitude using MSM
figure 
grid on
hold on
xlabel('Relative X distance [m]')
ylabel('Relative Y distance [m]')
zlabel('Relative Z distance [m]')
set(gca,'FontName','times')
makeSphsPicture_2craft( params.SPHS1, SPHS2_2, x0_H*1000,...
[0 0 0], [params.V1,params.V1]')
hold on
quiver3(0,0,0,L2_2(1,1)*nscale,L2_2(2,1)*nscale,L2_2(3,1)*nscale,'color',col_r,'Linewidth',3,'MaxHeadSize',1)
% title('Anti L Position')
c=colorbar;
c.Label.String = 'Surface Charge Density (nC/m^2)';
c.Label.FontSize = 14;
view(3)
axis tight
hold off
[ ~, ~, L1_2, L2_2, qs, overlapFlag] = ...
    multisphereFT( params.SPHS1, SPHS2_2, params.r, params.V,...
    params.C1, params.C2,params.COM1, params.COM2);

%% Computing Torques for all servicer Orientations (Stationary Target)

flag_solve4torques = false;

% Computing the Electrostatic Torques for different attitudes for SC 2
% This is computed using 3-2-1 Euler Angles with Yaw [-180 - 180], Pitch
% [-90 - 90], and Roll [-180 - 180] for n valus with in these ranges
% for n = 50, that is N = 125000
% if flag_solve4torques == true
%     n = 100;
%     data = All_Torques(params,n);
%     save('data')
% else
% %     load('data_15m_10V');
%     load('data')
% end
%%
% [data_anti,i,check] = find_anti_torque(L2_1,data);


% SPHS2_anti(1:3,:) = data_anti.C2*params.SPHS2(1:3,:);
% SPHS2_anti(4,:) = params.SPHS2(4,:);
% 
% % Plotting the Second attitude using MSM
% figure 
% grid on
% hold on
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, SPHS2_anti, x0_H*1000,...
% [0 0 0], params.V)
% hold on
% % title('Anti L Position')
% % c=colorbar;
% % c.Label.String = 'Surface Charge Density (nC/m^2)';
% % c.Label.FontSize = 14;
% view(3)
% axis tight
% hold off
% 
% [ ~, ~, L1_anti, L2_anti, qs, overlapFlag] = ...
%     multisphereFT( params.SPHS1, SPHS2_2, params.r, params.V,...
%     params.C1, params.C2,params.COM1, params.COM1);
% 



%%
%
data = 0;
[Xtot, Htot, Ttot, ttot, aterrtot, werrtot, ustot, Lrtot,L_e_tot,params,...
    normBN,reftot,rtot,flags] =...
    N_RW_sim(X0, I_RW, Iws, Gs0, N, K, P, tn, dt,params,data);

%
% % Show Plots
% flag_plots = true;
% % Make Video: takes a really long time, only use when necessary
% flag_vid = true;
% % Remember to change file name for each video
% filenam = 'Etorques_Tumble_testing_opposite_Torque';
%
flag_plots = true;
% Make Video: takes a really long time, only use when necessary
flag_vid = false;
% Remember to change file name for each video
filenam = 'E torques at various attitudes';

%%
if flag_plots == true
    ShowPlots(Xtot, Htot, Ttot, ttot, aterrtot, werrtot,ustot,Lrtot,L_e_tot,...
    params0,params,normBN,reftot,ttotal,x0_H,true) 
end
if flag_vid == true
    MakeVideo(Xtot,params,x0_H,filenam,rtot,Htot,L_e_tot)
end

% Used for Checking that the Electrostatic Torque computed at the 
% anti-torque location in the simulation is the same one as computed by the
% MSM

% L_check = L_e_tot(:,1);
% [data_anti,i_anti,tot] = find_anti_torque(L_check,data,params);
% check = dot(L_anti,L_check);


%% Major Functions

% Used for finding the attitude that results in a torque opposite that of H
function [data_anti,i_anti,anti_parr_check_min] = find_anti_torque(H,data)
%Must imput a non normalized H
% MAKE SURE ALL IN THE SAME FRAME
anti_parr_check_min = 1000;
data_anti = 0;
i_anti = 0;

for i = 1:length(data)
% If the spacecraft are overlapping, then set the value as NaN
    % if data{i}.overlapFlag == 0
        % If not, add the sum measure to tot
        L = data{i}.T_L;
        check = dot(L/norm(L),H/norm(H));

        % If the sum measure is closer to 0 than the previously
        % found minimum, set the minimum to be the current value.
        if check < anti_parr_check_min
            clear data_anti
            data_anti = data{i};
            i_anti = i;
            anti_parr_check_min = check;
        end
    % else
    %     tot(i,1:4) = NaN;
    % end
end

end




%%%%%%%%%%%%%%%%%%%Reaction Wheels Simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xtot, Htot, Ttot, ttot,aterrtot,werrtot,ustot,Lrtot,...
    L_e_tot,params,normBN,reftot,rtot,flags] = N_RW_sim(X0, I_RW, Iws, Gs0,...
    N, K, P, tn, dt,params,data)

X = X0; 
ttot = 0; 
Xtot = X0; Htot = []; Ttot = []; aterrtot = []; werrtot = []; ustot = [];
Lrtot = []; L_e_tot = []; normBN = []; reftot = []; rtot = []; flags = [];

r = params.r; V = params.V; C1 = params.C1; C2 = MRP2C(X(1:3));% C2 = params.C2;
r0 = r;
t = 0; 
refflag = 0; wheel_speed_flag = 0;
sig_RN = [0;0;0]; sig_anti = [ -1;0;0];
% sig_anti = [-0.0210321117924079,0.742001824556019,0.00229324345239622]';
Xprv = X0;

% [ ~, ~, ~, L20, ~, overlapFlag] = ...
%     multisphereFT( params.SPHS1, params.SPHS2, r, V, C1, C2,...
%     params.COM1, params.COM2);
% L = L20;
% L2 =  [ -0.000054008811829;0.003719000304351;-0.002708315879101];
while t < tn
%     if refflag == 0
%         sig_RN = [-1,0,0]';
%         refflag = 1;
%     elseif refflag == 1
%         if abs(X(1)-sig_RN(1)) < 0.01 
%             sig_RN = [0,-1,0]';
%             refflag = 2;
%         end
%     elseif refflag == 2
%         if abs(X(2)-sig_RN(2)) < 0.01 
%         sig_RN = [0,0,-1]';
%         refflag = 3;
%         end
%     end
 
    if refflag == 0
        i = 1;
        if wheel_speed_flag == 0
                X_sat = X0;
                wheel_speed_flag = 1;
        end
        while i < 4 
            if abs(X(i+6) - X_sat(i+6)) > 50
%             if mod(t,1000) == 0
%                 L = L2;
%                 BN = MRP2C(X(1:3));
%                 N_H = BN'*H;
%                 [data_anti,~,~] = find_anti_torque(N_H,data);
%                 sig_RN = DCM2MRP(data_anti.C2);
                sig_RN = sig_anti; 
                refflag = 1;
                sat_wheel = (i+6);
                X(1:3) = sig_RN;
                i = 4;
%                 X(1:3) = sig_anti;
%                 sig_RN = sig_anti;
%                 X(1:3) = [0,0,0]';
%                 sig_RN = [0,0,0]';
%                 r(1) = 10000;
                wheel_speed_flag = 0;
            end
            i = i + 1;
        end
    end
   
    if refflag == 1
        if abs(norm(X(1:3)) - norm(sig_RN)) < .001
            if wheel_speed_flag == 0
                X_sat = X;
                r(1) = r0(1);
                wheel_speed_flag = 1;
            end
            while i < 4   
                if abs(X_sat(sat_wheel) - X(sat_wheel)) > 50
%                if mod(t,1000) == 500
                    refflag = 2;
                    i = 4;
                    sig_RN = [0,0,0]'; 
                    X(1:3) = sig_RN;
%                     X(1:3) = [0,0,0]';
%                     sig_RN = [1,0,0]'; 
%                     X(1:3) = [1,0,0]';
                    r(1) = 10000;
                    wheel_speed_flag = 0;
                end
                
                i = i + 1;
            end
        end
    end
    
    if refflag == 2
        if abs(norm(X(1:3)) - norm(sig_RN)) < .001
            refflag = 0;
%             wheel_flag_speed = 0;
            r(1) = r0(1);
        end
    end
    


    reftot = [reftot,sig_RN];
    
    R_w_RN = [0,0,0]';
    R_w_RN_dot = [0,0,0]';
    
    sig_BN = X(1:3);
    w_BN = X(4:6);
    
    RN = MRP2C(sig_RN);
    BN = MRP2C(sig_BN);
    BR = BN*RN';
    sig_BR = C2MRP(BR);

    w_RN_B = -BR*R_w_RN;
    w_BR_B = w_BN - w_RN_B;
    w_dot_RN = BR*R_w_RN_dot;
    
    
    tt = 1;
    Om_mat = zeros(N,1);
    while tt < N+1
        Om_mat(tt) = X(6+tt);
        tt = tt + 1;
    end

    Gs = Gs0;
    hs = [];
    Hw = 0;
    for k = 1:N
    
        gsi = Gs(1:3,k);

        ws = dot(w_BN,gsi);

        hsi = Iws*(ws + Om_mat(k));
        hs = [hs;hsi];
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    
    end
    
%     Currently Assuming the Torques are given in the body frame of the
%     Target
    C2 = MRP2C(X(1:3));
    [ ~, ~, T_L1, T_L2, ~, overlapFlag] = ...
    multisphereFT( params.SPHS1, params.SPHS2, r, V, C1, C2,...
    params.COM1, params.COM2);
    
    B_L2 = C2*T_L2;
%         B_L2 = C2*T_L2;


% % 
%     if refflag == 1 %comment out
%         if wheel_speed_flag == 1
% %             B_L2 = L2; %comment out
%             B_L2 = C2*C1\L2;
%         end
% %     end
%     elseif wheel_speed_flag == 0 %comment out
% %         B_L2 = B_L2*0;%comment out
%     end %comment out
    

    L = B_L2;
        

    Lr = -K*sig_BR - P*(w_BR_B) + I_RW*((w_dot_RN) - tild(w_BN)*(w_RN_B)) + ...
    tild(w_BN)*(I_RW*w_BN + Gs*hs)-L;

    us = pinv(Gs)*-Lr;
    % Remove Control
%     us = us.*0;
    
    umax = .1;
    for i = 1:3
        if us(i) > umax
            us(i) = umax;
        elseif us(i) < -umax
            us(i) = -umax;
        end
    end

% % Setting maximum wheel speeds
%     Om_max = 300;
%     for i = 1:3
%         if X(6+i) > Om_max
%             if us(i) > 0
%                 us(i) = 0;
%             end
%         elseif X(6+i) < -Om_max
%             if us(i) < 0
%                 us(i) = 0;
%             end
%             
%         end
%     end


    aterr = norm(sig_BR);
    werr = norm(w_BR_B);

    
    Xdot = N_RW_EOM(X, I_RW, Iws, Gs0, N,us,L);

    k1 = Xdot*dt;
    Xdot = N_RW_EOM(X+k1/2, I_RW, Iws, Gs0, N,us,L);
    
    k2 = Xdot*dt;
    Xdot = N_RW_EOM(X+k2/2, I_RW, Iws, Gs0, N,us,L);

    k3 = Xdot*dt;
    Xdot = N_RW_EOM(X+k3, I_RW, Iws, Gs0, N,us,L);
    k4 = Xdot*dt;

    X = X + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set
    if norm(X(1:3))>1
       X(1:3) = -X(1:3)/norm(X(1:3))^2;
    end

    C2 = MRP2C(X(1:3));
    
    HB = I_RW*w_BN;
    
    H = HB + Hw;
    
 
    ttot = [ttot;t];
    t = t+dt;
    dX = X-Xprv;
    Xprv = X;
    Htot = [Htot,H];
    Xtot = [Xtot,X];
    aterrtot = [aterrtot,aterr];
    werrtot = [werrtot,werr];
    ustot = [ustot,us];
    Lrtot = [Lrtot,Lr];
%     L_e_tot = [L_e_tot,L];
    L_e_tot = [L_e_tot,T_L2];
    normBN = [normBN,norm(X(1:3))];
    rtot = [rtot,r];
    flags = [flags,[refflag;wheel_speed_flag]];


end
end

%Reaction Wheels Equations of Motion
function Xdot = N_RW_EOM(X, I_RW, Iws, Gs0, N,us,L)
sig_BN = X(1:3);
w_BN = X(4:6);

tt = 1;
Om_mat = zeros(N,1);
while tt < N+1
    Om_mat(tt) = X(6+tt);
    tt = tt + 1;
end

Gs = Gs0;
hs = [];
for k = 1:N
    
    gsi = Gs(1:3,k);  
    ws = dot(w_BN,gsi);  
    hsi = Iws*(ws + Om_mat(k));
    hs = [hs;hsi];

end

Bsig = ((1-norm(sig_BN)^2)*eye(3)+2*tild(sig_BN)+2*(sig_BN*sig_BN'));

sig_BN_dot = 1/4*Bsig*w_BN;

w_BN_dot = inv(I_RW)*(-tild(w_BN)*I_RW*w_BN - tild(w_BN)*(Gs*hs)-Gs*us+L);

Om_dot = us/Iws;
% Wheel speed saturation
% for i = 1:length(Om_mat)
%     if Om_mat(i) > 50
%         Om_dot(i) = 0;
%     end
% end

Xdot = [sig_BN_dot;w_BN_dot;Om_dot];
end



%%%%%%%%%%%%%%%Calculating Electrostatic Forces/Torques%%%%%%%%%%%%%%%%%%%%
function [ F1, F2, L1, L2, qs, overlapFlag] = ...
multisphereFT( SPHS1, SPHS2, r, V, C1, C2, COM1, COM2)
% Find force and torque between two MSM bodies
% Inputs:
% SPHS1 = [x1 x2 ...
%          y1 y2 ...
%          z1 z2 ...
%          R1 R2 ...]
%        (positions & radii of spheres composing body 1)
% SPHS2 = (same for body 2)
% C1/2 = DCM of body 1/2 rotation
% r = [xr;yr;zr] (position of body 2 origin relative to 1)
% V = [V1 V2] (voltages of bodies 1 and 2)
% COM1/2 = Center of mass to origin offset
%
% Outputs:
% F1, F2, L1, L2 = forces and torques (3x1 vectors) on body 1 and 2
% qs = vector of charges on each sphere

% assume zero rotations if not specified
if nargin < 5
    C1 = eye(3);
    C2 = eye(3);
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
elseif nargin < 6
    C2 = eye(3);
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
elseif nargin < 7
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
end

% Coulomb's constant [Nm^2/C^2]
k = constants.Kc;

% number of spheres in body 1, body 2, total
n1 = size(SPHS1,2);
n2 = size(SPHS2,2);
n = n1 + n2; 

% rotate positions of spheres according to DCMs
SPHS1t = SPHS1;
SPHS1t(1:3,:) = C1*(SPHS1(1:3,:)) + r*ones(1,n1); 
SPHS2(1:3,:) = C2*(SPHS2(1:3,:)); %

% build matrix with all spheres
SPHS = [SPHS1t SPHS2];

% Find Cinv matrix
% needs to evaluate every sphere to find each charge
sph1 = repmat(SPHS(1:3,:),1,size(SPHS,2));
sph2 = repelem(SPHS(1:3,:),1,size(SPHS,2));

relPos = sum((sph1-sph2).^2,1).^0.5; %vecnorm, but 1/10 the time. 
%Finds relative positions between each pair of spheres in the two bodies
%distance between each sphere: wil fill off-diagonal elements
relPos = reshape(relPos, n, n); 
Cinv = relPos + SPHS(4,:).*eye(n);

% check if any spheres are overlapping
lapCheck = relPos(n1+1:end, 1:n1) - 1.1*(SPHS1(4,:)+SPHS2(4,:)');
overlapFlag = 0;
if any(triu(lapCheck,1) < 0, 'all')
    overlapFlag = 1;
    %disp('Overlapping spheres')
%     figure; scatter3(SPHS1t(1,:),SPHS1t(2,:),SPHS1t(3,:));...
% hold on; scatter3(SPHS2(1,:),SPHS2(2,:),SPHS2(3,:),'filled');axis equal
% figure; makeSphsPicture_2craft( SPHS1t, SPHS2, ...
% [0 0 0], [0 0 0], V(1), eye(3), eye(3))
end


% figure; makeSphsPicture_2craft( SPHS1t, SPHS2,...
% [0 0 0], [0 0 0], V(1), eye(3), eye(3)); view(3)

Cinv = Cinv.^-1;

% Find charge on each sphere
Vs = [V(1)*ones(n1,1); V(2)*ones(n2,1)];
qs = (k*Cinv)\Vs;

% seperate charges into body 1 & 2
qs1 = qs(1:n1);
qs2 = qs(n1+1:end);

% Find force
% F = k*q1*q2*r/r^3

%compute distances 
sph1 = repmat(SPHS1t(1:3,:), 1, size(SPHS2,2));
sph2 = repelem(SPHS2(1:3,:), 1, size(SPHS1t,2));

relPos = sum((sph1-sph2).^2,1).^0.5;
relVec = sph1-sph2;

Qs1 = repmat(qs1',1,length(qs2));
Qs2 = repelem(qs2',1,length(qs1));

Feach = k*Qs1.*Qs2./relPos.^3;
Feach = Feach.*relVec; %add direction

F1 = sum(Feach,2); %matches to eps preision
F2 = -F1;

if norm(F1) >1
    %disp('High fores')
end
% Find torque
% L = r X f
%compute distances RELATIVE TO CENTER OF MASSES 
sph1 = repmat(SPHS1t(1:3,:) + COM1, 1, size(SPHS2,2));
sph22 = repelem(SPHS2(1:3,:) + COM2, 1, size(SPHS1t,2));

L1 = sum(cross(sph1,Feach),2); %matches L1 to eps
L2 = -sum(cross(sph22,Feach),2); % matches L2 to eps-ish

if norm(abs(L1)-abs(L2)) < eps() && nnz(V)>1
    %disp('Torques are equal')
end

end


function [Xtot, Htot, Ttot, ttot,aterrtot,werrtot,ustot,Lrtot,...
    L_e_tot,params,normBN,reftot,rtot,flags] = N_RW_sim_no_man(X0, I_RW, Iws, Gs0,...
    N, K, P, tn, dt,params,data,sig_RN)

X = X0; 
ttot = 0; 
Xtot = X0; Htot = []; Ttot = []; aterrtot = []; werrtot = []; ustot = [];
Lrtot = []; L_e_tot = []; normBN = []; reftot = []; rtot = []; flags = [];

r = params.r; V = params.V; C1 = params.C1; C2 = MRP2C(X(1:3));% C2 = params.C2;
r0 = 15;
t = 0; 
refflag = 0; wheel_speed_flag = 0;

Xprv = X0;

while t < tn
    


    reftot = [reftot,sig_RN];
    
    R_w_RN = [0,0,0]';
    R_w_RN_dot = [0,0,0]';
    
    sig_BN = X(1:3);
    w_BN = X(4:6);
    
    RN = MRP2C(sig_RN);
    BN = MRP2C(sig_BN);
    BR = BN*RN';
    sig_BR = C2MRP(BR);

    w_RN_B = -BR*R_w_RN;
    w_BR_B = w_BN - w_RN_B;
    w_dot_RN = BR*R_w_RN_dot;
    
    
    tt = 1;
    Om_mat = zeros(N,1);
    while tt < N+1
        Om_mat(tt) = X(6+tt);
        tt = tt + 1;
    end

    Gs = Gs0;
    hs = [];
    Hw = 0;
    for k = 1:N
    
        gsi = Gs(1:3,k);

        ws = dot(w_BN,gsi);

        hsi = Iws*(ws + Om_mat(k));
        hs = [hs;hsi];
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    
    end
    
    %Currently Assuming the Torques are given in the body frame of the
    %Target
%     [ ~, ~, T_L1, T_L2, ~, overlapFlag] = ...
%     multisphereFT( params.SPHS1, params.SPHS2, r, V, C1, C2,...
%     params.COM1, params.COM2);
%     
%     B_L2 = C2*T_L2;

%     B_L2 = [ -0.0045, -0.0084, 0.0009]';
    B_L2 = [0.00518*1.9;0.00386*.9;0.0000652/2];
%       B_L2 = [0.00518;0.00386;0.0000652];
%     

    L = B_L2;
        

    Lr = -K*sig_BR - P*(w_BR_B) + I_RW*((w_dot_RN) - tild(w_BN)*(w_RN_B)) + ...
    tild(w_BN)*(I_RW*w_BN + Gs*hs)-L;

    us = pinv(Gs)*-Lr;
    % Remove Control
    % us = us.*0;
    
    umax = .1;
    for i = 1:3
        if us(i) > umax
            us(i) = umax;
        elseif us(i) < -umax
            us(i) = -umax;
        end
    end

% % Setting maximum wheel speeds
%     Om_max = 50;
%     for i = 1:3
%         if X(6+i) > Om_max
%             if us(i) > 0
%                 us(i) = 0;
%             end
%         elseif X(6+i) < -Om_max
%             if us(i) < 0
%                 us(i) = 0;
%             end
%             
%         end
%     end


    aterr = norm(sig_BR);
    werr = norm(w_BR_B);

    
    Xdot = N_RW_EOM(X, I_RW, Iws, Gs0, N,us,L);

    k1 = Xdot*dt;
    Xdot = N_RW_EOM(X+k1/2, I_RW, Iws, Gs0, N,us,L);
    
    k2 = Xdot*dt;
    Xdot = N_RW_EOM(X+k2/2, I_RW, Iws, Gs0, N,us,L);

    k3 = Xdot*dt;
    Xdot = N_RW_EOM(X+k3, I_RW, Iws, Gs0, N,us,L);
    k4 = Xdot*dt;

    X = X + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set
    if norm(X(1:3))>1
       X(1:3) = -X(1:3)/norm(X(1:3))^2;
    end

    C2 = MRP2C(X(1:3));
    
    HB = I_RW*w_BN;
    
    H = HB + Hw;
    
 
    ttot = [ttot;t];
    t = t+dt;
    dX = X-Xprv;
    Xprv = X;
    Htot = [Htot,H];
    Xtot = [Xtot,X];
    aterrtot = [aterrtot,aterr];
    werrtot = [werrtot,werr];
    ustot = [ustot,us];
    Lrtot = [Lrtot,Lr];
    L_e_tot = [L_e_tot,L];
    normBN = [normBN,norm(X(1:3))];
    rtot = [rtot,r];
    flags = [flags,[refflag;wheel_speed_flag]];


end
end

function ShowPlots(Xtot, Htot, Ttot, ttot, aterrtot, werrtot,...
    ustot,Lrtot,L_e_tot,params0,params,normBN,reftot,ttotal,x0_H,flag) 
    col = [.466,.674,.188];
    if flag == true
        
        for i = 1:length(Htot)
            normH(i) = norm(Htot(:,i));
        end
        
        C2_intial = MRP2C(Xtot(1:3,1));
        
        intial_sph_loc = params.SPHS2*0;
        for i = 1:length(params.SPHS2)
            intial_sph_loc(1:3,i) = C2_intial*params.SPHS2(1:3,i);
            intial_sph_loc(4,i) = params.SPHS2(4,i);
        end
        C2_end = MRP2C(Xtot(1:3,end));
%         C2_end = MRP2C([-1;0;0]);
        final_sph_loc = params.SPHS2*0;
        for i = 1:length(params.SPHS2)
            final_sph_loc(1:3,i) = C2_end*params.SPHS2(1:3,i);
            final_sph_loc(4,i) = params.SPHS2(4,i);
        end

        figure
        plot(ttot,Xtot(1:3,:),'Linewidth',2)
        xlabel('Time (s)')
        ylabel('MRP Components')

        figure
        plot(ttot(2:end),normBN,'Linewidth',2)
        xlabel('Time (s)')
        ylabel('MRP Magnitude')

        figure
        subplot(2,2,1)
        plot(ttot,Xtot(1:3,:),'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('MRP Components','Fontsize',14)
        legend('sig\_1','sig\_2','sig\_3','Fontsize',14)
        xlim([0,ttotal])
        subplot(2,2,2)
        plot(ttot,Xtot(4:6,:),'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Angular Velocities (rad/s)','Fontsize',14)
        legend('w\_1','w\_2','w\_3','Fontsize',14)
        xlim([0,ttotal])
        subplot(2,2,3)
        plot(ttot,Xtot(7:9,:),'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Reaction Wheel Velocites (rad/s)','Fontsize',14)
        legend('Om\_1','Om\_2','Om\_3','Fontsize',14)
        subplot(2,2,4)
        xlim([0,ttotal])
        plot(ttot(2:end),ustot.*1000,'Linewidth',2)
        ylabel('Control Torque us (mNm)','Fontsize',14)
        xlabel('Time (s)','Fontsize',14)
        legend('us\_1','us\_2','us\_3','Fontsize',14)
        xlim([0,ttotal])


        figure
        plot(ttot,Xtot(7:9,:),'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Reaction Wheel Velocites (rad/s)','Fontsize',14)
        legend('Om\_1','Om\_2','Om\_3','Fontsize',14)
        
        figure
        plot(ttot(2:end),reftot(1:3,:),'Linewidth',2)
        xlabel('Time (s)')
        ylabel('Desired Attitude (MRPs)')

        figure
        subplot(2,1,1)
        semilogy(ttot(2:end),aterrtot,'Linewidth',2)
        ylabel('Magnitude sigma B/R ','Fontsize',14)
        subplot(2,1,2)
        semilogy(ttot(2:end),werrtot,'Linewidth',2)
        xlabel('Time (s)','Fontsize',14)
        ylabel('Magnitude omega B/R ','Fontsize',14)
        
        figure
%         subplot(2,1,1)
%         plot(ttot(2:end),Htot)
%         ylabel('Angular Momentum')
%         xlabel('Time (s)')  
        plot(ttot/3600,Xtot(7,:)/2/pi*60,'Linewidth',2)
        hold on
        plot(ttot/3600,Xtot(8,:)/2/pi*60,'Linewidth',2)
        plot(ttot/3600,Xtot(9,:)/2/pi*60,'color',col,'Linewidth',2)
        yline(4000,'--')
        xlabel('Time (hr)','Fontsize',14,'Fontname','Times New Roman')
        ylabel('RW Velocity (rpm)','Fontsize',14,'Fontname','Times New Roman')
        legend('X','Y','Z','Fontsize',12,'Fontname','Times New Roman','Location','northwest')
%         subplot(2,1,2)
        figure
        plot(ttot(2:end)/3600,Htot(1,:),'Linewidth',2)
        hold on
        plot(ttot(2:end)/3600,Htot(2,:),'Linewidth',2)
        plot(ttot(2:end)/3600,Htot(3,:),'color',col,'Linewidth',2)
        legend('X','Y','Z','Fontsize',14,'location','northwest')
%         xline(38120,'--')
        ylabel('Angular Momentum (Nms)','Fontsize',14,'Fontname','Times New Roman')
        xlabel('Time (hr)','Fontsize',14,'Fontname','Times New Roman')
%         ylim([0,350])
        xlim([0,ttot(end)/3600])
        box off
        
        figure
        plot(ttot(2:end)/3600,normH,'Linewidth',2)
        hold on
%         xline(38120,'--')
        ylabel('Angular Momentum (Nms)','Fontsize',14,'Fontname','Times New Roman')
        xlabel('Time (hr)','Fontsize',14,'Fontname','Times New Roman')
%         ylim([0,350])
        xlim([0,ttot(end)/3600])
        box off

        figure
        plot(ttot(2:end),L_e_tot,'Linewidth',2)
        ylabel('Electrostatic Torques (Nm)','Fontsize',14)
        xlabel('Time (s)','Fontsize',14)
        legend('L\_x','L\_y','L\_z','Fontsize',14)

        %Initial position plot
        figure 
        grid on
        hold on
        xlabel('Relative X distance [m]')
        ylabel('Relative Y distance [m]')
        zlabel('Relative Z distance [m]')
        set(gca,'FontName','times')
        makeSphsPicture_2craft( params.SPHS1, intial_sph_loc, x0_H*1000,...
        [0 0 0], params.V)
%         quiver3(0,0,0,40*Htot(1,1)/100,40*Htot(2,1)/100,40*Htot(3,1)/100,'r','Linewidth',3)
        quiver3(0,0,0,L_e_tot(1,1)*10000,L_e_tot(2,1)*10000,L_e_tot(3,1)*10000,'k','Linewidth',3)
%         quiver3(0,0,0,10e2*L_e_tot(1,1),10e2*L_e_tot(2,1),15e2*L_e_tot(3,1),'b','Linewidth',3)
        axis tight
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
%         c.Label.FontSize = 14;
%         title('Initial Positions')
        view(3)
        
        %Final position plot
        figure 
        grid on
        hold on
        xlabel('Relative X distance [m]')
        ylabel('Relative Y distance [m]')
        zlabel('Relative Z distance [m]')
        set(gca,'FontName','times')
        makeSphsPicture_2craft( params.SPHS1, final_sph_loc, x0_H*1000,...
        [0 0 0], params.V)
%         quiver3(0,0,0,Htot(1,end)/5,Htot(2,end)/5,Htot(3,end)/5,'k','Linewidth',3)
        quiver3(0,0,0,L_e_tot(1,end)*10000,L_e_tot(2,end)*10000,L_e_tot(3,end)*10000,'k','Linewidth',3)
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
%         c.Label.FontSize = 14;
        axis tight
%         title('Final Positions')
        view(3)
    end

end


function MakeVideo(Xtot,params,x0_H,filenam,rtot,Htot,L_e_tot)
n = (length(Xtot)-1)/1;
inter_sph_loc = params.SPHS2*0;
writerObj = VideoWriter(filenam,'MPEG-4');
v.Quality = 100;
open(writerObj);
heights = zeros(n,1);
lengths = zeros(1,n);
for j = 2:n
    C2 = MRP2C(Xtot(1:3,j));
    
    for i = 1:length(params.SPHS2)
        inter_sph_loc(1:3,i) = C2*params.SPHS2(1:3,i);
        inter_sph_loc(4,i) = params.SPHS2(4,i);
    end
 
    f = figure(20);
%     f.Position = [3000 800 750 750];
    makeSphsPicture_2craft(params.SPHS1, inter_sph_loc,...
    rtot(:,j), [0 0 0], params.V); %x0_H*1000
%     quiver3(0,0,0,Htot(1,j)/5,Htot(2,j)/5,Htot(3,j)/5,'k','Linewidth',3)
    quiver3(0,0,0,L_e_tot(1,j)*10000,L_e_tot(2,j)*10000,L_e_tot(3,j)*10000,'k','Linewidth',3)

    view(3)
    axis tight
%     drawnow
%     ax = gca;
%     ax.Units = 'pixels';
    F(j-1) = getframe(f);
%     F(j-1) = getframe;
    size_F = size(F(j-1).cdata);
    heights(j-1) = size_F(1);
    lengths(j-1) = size_F(2);

    if j ~= n
        clf(figure(20))
    end
end
max_height = max(heights);
max_length = max(lengths);

F = F(2:end);
for k = 1:length(F)
    if lengths(k) < max_length
        F(k).cdata(:,lengths(k)+1:max_length,:) = 0;
    end
    if heights(k) < max_height
        F(k).cdata(heights(k)+1:max_height,:,:) = 0;
    end
    writeVideo(writerObj,F(k))
end
close(writerObj);
end

