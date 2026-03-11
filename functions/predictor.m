function [metric, omegaTargetF,params] = predictor(polpar,t0,inicond,refTraj,params)
% Based on the alpha and beta parameters tested by the solver, this
% function predicts the final rotation rate of the target frame

%% Extract polynomial params
%polpar = [1   polpar(1)  1 polpar(2)];
% alpha = polpar(1:2);
% beta  = polpar(3:4);

polpar = [polpar(1)  polpar(2)];
alpha  = polpar(1);
beta   = polpar(2);

%% Extract initial position and velocity
DCM_TH         = quat2dcm(inicond(20:23)');            % Target object DCM
rho_ST_T       = DCM_TH*(inicond(1:3)-inicond(14:16)) - params.trajOffset; % position servicer to the target
rhoDot_ST_T    = DCM_TH*(inicond(4:6)-inicond(17:19));
X0            = rho_ST_T(1)*1000;
Y0            = rho_ST_T(2)*1000;
Z0            = rho_ST_T(3)*1000;
% VX0           = rhoDot_ST_T(1)*1000;
% VY0           = rhoDot_ST_T(2)*1000;
% VZ0           = rhoDot_ST_T(3)*1000;
VX0           = params.refTrajInt.dX(t0);
VY0           = params.refTrajInt.dY(t0);
VZ0           = params.refTrajInt.dZ(t0);

%% Compute new trajectory in T frame
refTraj = polynomialinterpolant(X0,Y0,Z0,VX0,VY0,VZ0,alpha,beta,refTraj);

% Store in params structure
params.refTraj       = refTraj.refTraj'; %TRAJECTORY IN T FRAME
params.tRef          = refTraj.tRef;     %timestamps for the reference trajectory...

% Express as griddedInterpolants for sweet performance
params.refTrajInt.X  = griddedInterpolant(params.tRef, params.refTraj(:,1) + params.trajOffset(1)*1000);
params.refTrajInt.Y  = griddedInterpolant(params.tRef, params.refTraj(:,2) + params.trajOffset(2)*1000);
params.refTrajInt.Z  = griddedInterpolant(params.tRef, params.refTraj(:,3) + params.trajOffset(3)*1000);
params.refTrajInt.dX = griddedInterpolant(params.tRef, params.refTraj(:,4));
params.refTrajInt.dY = griddedInterpolant(params.tRef, params.refTraj(:,5));
params.refTrajInt.dZ = griddedInterpolant(params.tRef, params.refTraj(:,6));

%% Integrate equations in HCW
% numerical integration of Cartesian State Vector
% STRUCTURE FOR SV:
% [1:3]-position 1          |   [14:16]-position 2
% [4:6]-velocity 1          |   [17:19]-Velocity 2
% [7:10]-quat 1             |   [20:23]-quat 2
% [11:13]-omega1            |   [24:26]-omega 2

% Parameters
tSpan          = [t0,params.tRef(end)];
options        = odeset('Maxstep', 100, 'RelTol', 1e-8, 'AbsTol', 1e-8);
params.perturb = 'electro';
params.srpPert = 'on';
params.sHat0   = [1 0 0];
params.sHat0   = params.sHat0/norm(params.sHat0);
params.control     = 'on';
params.feedForward = 'on'; %off
params.navNoise    = 'off';

% Integrate
[~, statesPerturbed] = ode15s(@(t,y) hcwDynamics6DOF_MSM_C2(t,y,params),tSpan,inicond,options);

% Extract angular velocity
omegaTarget = statesPerturbed(:, 24:26);
omegaTargetF= vecnorm(omegaTarget(end,:)*180/pi);

% Check visual cone condition
anglemax = max(atan(sqrt(params.refTraj(:,2).^2 + params.refTraj(:,3).^2)./params.refTraj(:,1))*180/pi);
if anglemax > params.conemax
    metric = omegaTargetF * (1+10*(anglemax-params.conemax));
else 
    metric = omegaTargetF;
end

%% Debug
fprintf('Final target rot rate (electrostatic perturbed): %f deg/s (maxang = %f deg) \n', omegaTargetF, anglemax)

% pos1Pert = statesPerturbed(:,1:3)*1000;
% pos2Pert = statesPerturbed(:,14:16)*1000;
% 
% figure %Relative position plot
% plot3(pos1Pert(:,1),pos1Pert(:,2),pos1Pert(:,3),'linewidth',2)
% grid on
% hold on
% plot3(pos2Pert(:,1),pos2Pert(:,2),pos2Pert(:,3),'linewidth',2)
% xlabel('Relative X distance [m]')
% ylabel('Relative Y distance [m]')
% zlabel('Relative Z distance [m]')
% scatter3(pos1Pert(1,1),pos1Pert(1,2),pos1Pert(1,3),'o','filled')
% legend('Servicer, perturbed','Target, perturbed','AutoUpdate','off')
% set(gca,'FontName','times')
% makeSphsPicture_2craft( params.SPHS1, params.SPHS2, pos1Pert(end,:), pos2Pert(end,:), params.V(1), quat2dcm(statesPerturbed(end,7:10))', quat2dcm(statesPerturbed(end,20:23))')
