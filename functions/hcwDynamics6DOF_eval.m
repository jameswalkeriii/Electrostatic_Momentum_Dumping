function [C_H, deltaRho_T, deltaRhoPrime_T, aElec_H, F1_H, F2_H, L1_H, L2_H, rho_ST_T, rhoDot_ST_T, Qs] = hcwDynamics6DOF_eval(t, u, params)
% Just returns control-relevant parameters for evaluation

% STRUCTURE FOR U: 
% [1:3]-position 1          |   [14:16]-position 2
% [4:6]-velocity 1          |   [17:19]-Velocity 2
% [7:10]-quat 1             |   [20:23]-quat 2
% [11:13]-omega1            |   [24:26]-omega 2

n = params.n;

if strcmp(params.control, 'on')

    desPos_T = [params.refTrajInt.X(t); params.refTrajInt.Y(t); params.refTrajInt.Z(t)];
    desVel_T = [params.refTrajInt.dX(t); params.refTrajInt.dY(t); params.refTrajInt.dZ(t)];
    
    
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
    
    desStates_T = [desPos_T; desVel_T]/1000;  %find reference trajectory [km]
    
    DCM_TH = quat2dcm(u(20:23)'); %Target object DCM

    rho_ST_T = DCM_TH*(u(1:3)-u(14:16)); %position servicer to the target
    rhoDot_ST_T = DCM_TH*(u(4:6)-u(17:19));

    deltaRho_T = ( rho_ST_T - desStates_T(1:3)); % position error from current relative pos to desired
    deltaRhoPrime_T = ( rhoDot_ST_T - desStates_T(4:6));
    
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
        disp('Maxed out C on eval')
    end
else
    C_H = zeros(1,3); %for case with no control
end

if strcmp(params.perturb , 'electro') %compute electrostatic forces/torques
    DCM_SH = quat2dcm(u(7:10)');
    DCM_TH = quat2dcm(u(20:23)');
    r_ST_H = (u(1:3)-u(14:16))*1000; % distance from servicer to target in Hill [m]
    V1     = params.V1(t);
    V2     = params.V2(t);
    [ F1_H, F2_H, L1_H, L2_H, Qs] = multisphereFT( params.SPHS1, params.SPHS2, r_ST_H, [V1;V2], DCM_SH', DCM_TH', params.COM1, params.COM2 ); %computes forces and torques

    if norm(F1_H) > 1 || norm(F2_H) > 1  %check for higher forces than expected
        F1_H = F1_H/norm(F1_H);
        F2_H = F2_H/norm(F2_H);
        fprintf('Forces greater than threshold, likely singular case at t= %d \n',t)
        disp(F1_H)
    end
    
    aElec_H(1:3) = F1_H/(params.m1*1000);% Acceleration due to electrostatic perturbations in km/s
    aElec_H(4:6) = F2_H/(params.m2*1000); 
    LElec_H(1:6) = [L1_H; L2_H]; %Add torque
else 
    aElec_H = zeros(1,6);
    F1_H = zeros(1,3);
    F2_H = zeros(1,3);
    L1_H = zeros(1,3); 
    L2_H = zeros(1,3);
    LElec_H(1:6) = zeros(1,6); %Added to torque computation later
end

aElec_H = aElec_H';

end

