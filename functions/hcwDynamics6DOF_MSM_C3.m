function dy = hcwDynamics6DOF_MSM_C3(t, u, params)
% STRUCTURE FOR U: 
% [1:3]-position 1          |   [14:16]-position 2
% [4:6]-velocity 1          |   [17:19]-Velocity 2
% [7:10]-quat 1             |   [20:23]-quat 2
% [11:13]-omega1            |   [24:26]-omega 2
% 1 IS SERVICER                 2 IS TARGET
n = params.n;
if strcmp(params.control, 'on')
    desPos_T = [params.refTrajInt.X(t); params.refTrajInt.Y(t); params.refTrajInt.Z(t)];
    desVel_T = [params.refTrajInt.dX(t); params.refTrajInt.dY(t); params.refTrajInt.dZ(t)];
    desStates_T = [desPos_T; desVel_T]/1000;  % find reference trajectory [km]
    DCM_TH = quat2dcm(u(20:23)'); % Target object DCM
    rho_ST_T = DCM_TH*(u(1:3)-u(14:16)); % position servicer to the target        
    rhoDot_ST_T = DCM_TH*(u(4:6)-u(17:19));
    if strcmp(params.navNoise, 'on')
        rho_ST_T = rho_ST_T + [params.navNoiseVal.X(t); params.navNoiseVal.Y(t); params.navNoiseVal.Z(t)]; % assume noise on the order of 10 cm
        rhoDot_ST_T = rhoDot_ST_T + [params.navNoiseVal.dX(t); params.navNoiseVal.dY(t); params.navNoiseVal.dZ(t)]; %velocity error of ~1 cm/s
    end
    deltaRho_T = 1000*(rho_ST_T - desStates_T(1:3)); % position error from current relative pos to desired
    deltaRhoPrime_T = 1000*(rhoDot_ST_T - desStates_T(4:6));
    deltaAccel_H = zeros(3,1);
    % Feedforward control using target body rotation rate/ acceleration
    % between reference point in target frame and Hill frame
    if strcmp(params.feedForward, 'on')
        desPos_H = DCM_TH*desStates_T(1:3);
        desVel_H = DCM_TH*desStates_T(4:6);
        % Feed forward on gravitational acceleration delta between
        % reference position and servicer
        deltaAccel_H(1) = (3*n^2*u(1) + 2*n*u(5)) - (3*n^2*desPos_H(1) + 2*n*desVel_H(2));
        deltaAccel_H(2) = (-2*n*u(4)) - (-2*n*desVel_H(1));
        deltaAccel_H(3) = (-n^2*u(3)) - (-n^2*desPos_H(3));
        % Target body rotation-induced acceleration
        rDes_TS_H = DCM_TH'*desPos_T;
        omega_TH = u(24:26);
        rDotDes_TS_H = DCM_TH'*desVel_T;
        rhoDot_TS_H = rDotDes_TS_H + cross(omega_TH, rDes_TS_H);
        rhoDDot_TS_H = cross(omega_TH, rhoDot_TS_H);
        deltaAccel_H = deltaAccel_H + rhoDDot_TS_H;
    end
    C_H = -DCM_TH'*(params.K1*deltaRho_T + params.K2*deltaRhoPrime_T) + deltaAccel_H;   % Control effort in H frame coordinates
    C_H = C_H/1000;  % put in km/s^2
    if norm(C_H) > params.maxC  % controller saturation check
        C_H = params.maxC*C_H/norm(C_H);
        %disp('Maxed out C')
    end
else
    C_H = zeros(1,3); %for case with no control
end
srpForce = [0 0 0];
srpTorque = [0 0 0];
if strcmp(params.srpPert, 'on')
    % compute sHat; use mean orbital motion (n) and time and rotx()
    sHat = params.sHat0 * roty(n*t*180/pi);
    sHat = sHat/norm(sHat);
    srPressure = 2.7e-6; % N/m https://www.hindawi.com/journals/ijae/2015/928206/tab2/
    for jj = 1:length(params.A)
        if dot(sHat, params.nHat(jj,:)) > 0
            cosTheta = dot(sHat, params.nHat(jj,:))/(norm(sHat)*norm(params.nHat(jj,:)));
            force = -srPressure*((1-params.spec(jj))*sHat + 2*(params.spec(jj)*cosTheta + 1/3*params.diffuse(jj))*params.nHat(jj,:))*cosTheta*params.A(jj);
            srpForce = srpForce + force;
            srpTorque = srpTorque + cross(params.C(jj,:), force); 
        end
    end
end
srpAccel = srpForce/(params.m2*1000); % in km/s^2
if strcmp(params.perturb, 'electro') %compute electrostatic forces/torques
    DCM_SH = quat2dcm(u(7:10)');
    DCM_TH = quat2dcm(u(20:23)');
    r_ST_H = (u(1:3)-u(14:16))*1000; % distance from servicer to target in Hill [m]
    V1     = params.V1(t);
    V2     = params.V2(t);
    [ F1_H, F2_H, L1_H, L2_H] = multisphereFT( params.SPHS1, params.SPHS2, r_ST_H, [V1;V2], DCM_SH', DCM_TH', params.COM1, params.COM2 ); %computes forces and torques
    if norm(F1_H) > 1 || norm(F2_H) > 1  %check for higher forces than expected
        F1_H = F1_H/norm(F1_H);
        F2_H = F2_H/norm(F2_H);
        fprintf('Forces greater than threshold, likely singular case at t= %d \n',t)
        disp(F1_H)
    end
    aElec_H(1:3) = F1_H/(params.m1*1000);% Acceleration due to electrostatic perturbations in km/s
    aElec_H(4:6) = F2_H/(params.m2*1000); 
    LElec_H(1:6) = [L1_H; L2_H];%Add torque
else 
    aElec_H = zeros(1,6);
    LElec_H(1:6) = zeros(1,6); %Added to torque computation later
end
%% SC1 translation-SERVICER
x = u(1); 
% y = u(2);
z = u(3);
xdot = u(4); 
ydot = u(5); 
% zdot = u(6);
xddot = 3*n^2*x+2*n*ydot + aElec_H(1) + C_H(1);
yddot = -2*n*xdot + aElec_H(2) + C_H(2);
zddot = -n^2*z + aElec_H(3) + C_H(3);
dx_S = [u(4);            u(5);           u(6);  xddot;  yddot; zddot ];
%% Find Target translation
x = u(14); 
% y = u(15);
z = u(16);
xdot = u(17); ydot = u(18); zdot = u(19);
xddot = 3*n^2*x+2*n*ydot + aElec_H(4) + srpAccel(1);
yddot = -2*n*xdot + aElec_H(5) + srpAccel(2);
zddot = -n^2*z + aElec_H(6) + srpAccel(3);
dx_T = [xdot; ydot; zdot; xddot; yddot; zddot];
%% rotational EOMs
beta_SH = u(7:10); %EP of Servicer
beta_SH = beta_SH/norm(beta_SH);
beta_TH = u(20:23);
beta_TH = beta_TH/norm(beta_TH); %quaterions of SC2
% omega1 = u(11:13); %rotational rate of SC1
omega_T = u(24:26); %rotational rate of target
B =@(beta) [-beta(2) -beta(3) -beta(4);
            beta(1) -beta(4) beta(3);
            beta(4) beta(1) -beta(2);
            -beta(3) beta(2) beta(1)];
betaDot2 = .5*B(beta_TH)*omega_T;
omegaTDot = params.I2\(cross(-omega_T, params.I2*omega_T) + LElec_H(4:6)' + srpTorque') ;
if strcmp(params.trackingOmega,'locked')
    omega_S = omega_T; %assuming att alg tracks target's rotational rates
    omegaSDot = omegaTDot;
elseif strcmp(params.trackingOmega,'follow')
    omega_S = u(11:13);
    delBeta = quatDiff(beta_SH, beta_TH);
    delOmega = omega_S - omega_T;
    controlTorque = -params.K_att*delBeta(2:4) - params.P_att*delOmega + params.I1*(omegaTDot-cross(omega_S, omega_T)) + tilde(omega_S)*params.I1*omega_S;
    omegaSDot = params.I1\(cross(-omega_S, params.I1*omega_S) + LElec_H(1:3)' + controlTorque);
elseif strcmp(params.trackingOmega,'refTrack')
    omega_S = u(11:13);
    DCM_SH = quat2dcm(u(7:10)');
    r_ST_H = (u(1:3)-u(14:16))*1000; % distance from servicer to target in Hill [m]
    guidedAngle = rotx(params.refTrajInt.e1(t)); % angle picked by the guidance algo
    LOSvec = r_ST_H; % line of sight from servicer to target
    dockingPortVec = DCM_SH*[1 0 0]'; % orientation of docking port on servicer (servicer attitude, basically)
    constrainedAngle = vec2dcm(LOSvec, dockingPortVec);
    beta_des = dcm2quat(constrainedAngle*guidedAngle); % desired quaternion
    delBeta = quatDiff(beta_SH, beta_des); % error in quaternion
%     omega_des = [params.refTrajInt.de1(t) + omega_T(1); % desired rotational rates
%                  omega_T(2); 
%                  omega_T(3)];
%              
    delOmega = omega_S - omega_T;%omega_des;
    controlTorque = -params.K_att*delBeta(2:4) - params.P_att*delOmega + params.I1*(omegaTDot-cross(omega_S, omega_T)) + tilde(omega_S)*params.I1*omega_S;
%     servicerSphs = DCM_SH*params.SPHS1(1:3,:) + r_ST_H;
%     targetSphs = DCM_TH*params.SPHS2(1:3,:);
%     figure
%     scatter3(servicerSphs(1,:), servicerSphs(2,:), servicerSphs(3,:), 'o')
%     hold on
%     scatter3(targetSphs(1,:), targetSphs(2,:), targetSphs(3,:), 'o', 'filled')
%     beta_des = dcm2quat(constrainedAngle*guidedAngle'); % desired quaternion
% 
%     servicerSphs = DCM_SH*params.SPHS1(1:3,:);
%     desServeSphs = quat2dcm(beta_des)*params.SPHS1(1:3,:);
%     figure
%     scatter3(servicerSphs(1,:), servicerSphs(2,:), servicerSphs(3,:), 'o')
%     hold on
%     scatter3(desServeSphs(1,:), desServeSphs(2,:), desServeSphs(3,:), 'o', 'filled')
%     axis equal
%     controlTorque = -params.K_att*delBeta(2:4) - params.P_att*delOmega + params.I1*(omegaTDot-cross(omega_S, omega_des)) + tilde(omega_S)*params.I1*omega_S;
    omegaSDot = params.I1\(cross(-omega_S, params.I1*omega_S) + LElec_H(1:3)' + controlTorque);
end
if norm(omegaSDot)>0.5
    disp('too much spinny')
end
betaDot1 = .5*B(beta_SH)*omega_S;
dy = [dx_S; betaDot1; omegaSDot; dx_T; betaDot2; omegaTDot];
end