function dy = rdzDynGuidCtrl_r2(t, u, params)
% STRUCTURE FOR U: 
% [1:3]-position 1          |   [14:16]-position 2
% [4:6]-velocity 1          |   [17:19]-Velocity 2
% [7:10]-quat 1             |   [20:23]-quat 2
% [11:13]-omega1            |   [24:26]-omega 2
% 1 IS SERVICER                 2 IS TARGET
n = params.n;
sigma = DCM2MRP(eye(3));
% SRP contributions
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
            forceVec = -srPressure*((1-params.spec(jj))*sHat + 2*(params.spec(jj)*cosTheta + 1/3*params.diffuse(jj))*params.nHat(jj,:))*cosTheta*params.A(jj);
            srpForce = srpForce + forceVec;
            srpTorque = srpTorque + cross(params.C(jj,:), forceVec); 
        end
    end
end
if strcmp(params.control, 'on')
    desPos_T = [params.refTrajInt.X(t); params.refTrajInt.Y(t); params.refTrajInt.Z(t)]; % in [m]
    desVel_T = [params.refTrajInt.dX(t); params.refTrajInt.dY(t); params.refTrajInt.dZ(t)];
%     if strcmp(params.guide, 'on')
%         DCM_SH = quat2dcm(u(7:10)');
%         DCM_TH = quat2dcm(u(20:23)');
%         r = params.searchPositions * norm(desPos_T);
%         L2 = zeros(size(r));
%         targetSphs = params.SPHS2;
%         targetSphs(1:3,:) = DCM_TH * params.SPHS2(1:3,:);
%         serveSphs = params.SPHS1;
%         
%         for ii = 1:length(r)
%             rtest = r(:,ii)';% * C1';
%             C2 = vec2dcm(DCM_TH*[1 0 0]', rtest');
%             serveSphs(1:3,:) = C2' *  DCM_SH * rotx(0) *  params.SPHS1(1:3,:); %Orients Servicer along right axis to point at docking point
%             [ ~, ~, ~, L2(:,ii), ~, ~] = multisphereFT( serveSphs, targetSphs,  rtest', params.V, eye(3), eye(3), params.COM1, params.COM2 );
%         end
%        
%         LNormNet = vecnorm(L2 + repmat(srpTorque, size(L2,2),1)' - repmat(params.Ldes, size(L2,2),1)'); %norm of net torque on target
%         [~, idx] = min(abs(LNormNet)); %find point that comes closest to meeting desired torque
%         
% %         % UPDATE SERVICER ATTITUDE
% %         % find MRP
% %         sigma = DCM2MRP(C1');
%         
%         % define new desStates_T (position and velocity)
%         desPos_T = r(:,idx);
%         desVel_T = DCM_SH*desVel_T;
%         
%     end
    desStates_T = [desPos_T; desVel_T]/1000;  %find reference trajectory [km]
    DCM_TH = quat2dcm(u(20:23)'); %Target object DCM
    rho_ST_T = DCM_TH*(u(1:3)-u(14:16)); %position servicer to the target        
    rhoDot_ST_T = DCM_TH*(u(4:6)-u(17:19));
    if strcmp(params.navNoise, 'on')
        rho_ST_T = rho_ST_T + [params.navNoiseVal.X(t); params.navNoiseVal.Y(t); params.navNoiseVal.Z(t)]; % assume noise on the order of 10 cm
        rhoDot_ST_T = rhoDot_ST_T + [params.navNoiseVal.dX(t); params.navNoiseVal.dY(t); params.navNoiseVal.dZ(t)]; %velocity error of ~1 cm/s
    end
    deltaRho_T = 1000*(rho_ST_T - desStates_T(1:3)); % position error from current relative pos to desired
    deltaRhoPrime_T = (rhoDot_ST_T - desStates_T(4:6));
    deltaAccel_H = zeros(3,1);
%    Feedforward terms for relative gravitational acceleration
    if strcmp(params.feedForward, 'on')
        % Find keplerian relative accelerations
        desPos_H = DCM_TH*desStates_T(1:3);
        desVel_H = DCM_TH*desStates_T(4:6);
        deltaAccel_H(1) = (3*n^2*u(1) + 2*n*u(5)) - (3*n^2*desPos_H(1) + 2*n*desVel_H(2));
        deltaAccel_H(2) = (-2*n*u(4)) - (-2*n*desVel_H(1));
        deltaAccel_H(3) = (-n^2*u(3)) - (-n^2*desPos_H(3));
        % Track general target body rotation:
        rDes_TS_H = DCM_TH'*desPos_T;
        omega_TH = u(24:26);
        rDotDes_TS_H = DCM_TH'*desVel_T;
        rhoDot_TS_H = rDotDes_TS_H + cross(omega_TH, rDes_TS_H);
        rhoDDot_TS_H = cross(omega_TH, rhoDot_TS_H); % diff(desVel_T = [params.refTrajInt.dX(t); params.refTrajInt.dY(t); params.refTrajInt.dZ(t)];)
        deltaAccel_H = deltaAccel_H + rhoDDot_TS_H;
    end
    C_H = -DCM_TH'*(params.K1*deltaRho_T + params.K2*deltaRhoPrime_T) + deltaAccel_H;   % Control effort in H frame coordinates
    C_H = C_H/1000;  % put in km/s^2
    if norm(C_H) > params.maxC %controller saturation
        C_H = params.maxC*C_H/norm(C_H);
        %disp('Maxed out C')
    end
else
    C_H = zeros(1,3); %for case with no control
end
srpAccel = srpForce/(params.m2*1000); % in km/s^2
if strcmp(params.perturb, 'electro') %compute electrostatic forces/torques
    DCM_SH = quat2dcm(u(7:10)');
    DCM_TH = quat2dcm(u(20:23)');
    r_ST_H = (u(1:3)-u(14:16))*1000; % position of target relative to servicer in Hill [m]
    V1     = params.V1(t);
    V2     = params.V2(t);
    [ F1_H, F2_H, L1_H, L2_H] = multisphereFT( params.SPHS1, params.SPHS2, r_ST_H, [V1;V2], DCM_SH', DCM_TH', params.COM1, params.COM2 ); %computes forces and torques
    if norm(F1_H) > 1 || norm(F2_H) > 1  %check for higher forces than expected
        F1_H = F1_H/norm(F1_H);
        F2_H = F2_H/norm(F2_H);
        fprintf('Forces greater than threshold (%.2f), likely singular case at t= %.3d \n', norm(F2_H), t)
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
z = u(3);
xdot = u(4); 
ydot = u(5); 
xddot = 3*n^2*x+2*n*ydot + aElec_H(1) + C_H(1);
yddot = -2*n*xdot + aElec_H(2) + C_H(2);
zddot = -n^2*z + aElec_H(3) + C_H(3);
dx_S = [u(4);            u(5);           u(6);  xddot;  yddot; zddot ];
%% Find Target Hill-frame translation
x = u(14); z = u(16);
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
omega_S = u(11:13); %rotational rate of SC1
omega_T = u(24:26); %rotational rate of target
B =@(beta) [-beta(2) -beta(3) -beta(4);
            beta(1) -beta(4) beta(3);
            beta(4) beta(1) -beta(2);
            -beta(3) beta(2) beta(1)];
betaDot2 = .5*B(beta_TH)*omega_T;
omegaTDot = params.I2\(cross(-omega_T, params.I2*omega_T) + LElec_H(4:6)' + srpTorque') ;
if strcmp(params.attControl,'off')
    omega_S = omega_T; %assuming att alg tracks target's rotational rates
    omegaSDot = omegaTDot;
elseif strcmp(params.attControl,'on')% align servicer to face docking port on target
    % use MRP feedback
    Lcontrol = -.02*sigma - .03*eye(3)* omega_S - LElec_H(1:3)';
    omega_S = u(11:13);
    omegaSDot = params.I1\(cross(-omega_S, params.I1*omega_S) + LElec_H(1:3)' + Lcontrol);
end
betaDot1 = .5*B(beta_SH)*omega_S;
dy = [dx_S; betaDot1; omegaSDot; dx_T; betaDot2; omegaTDot];
end
