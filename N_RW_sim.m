function results = N_RW_sim(X0_serv,X0_deb, Iws, Gs0, K, P, Ttot,params,data,K_h)

I_RW = params.servicer.S_MI; % The moment of inertia of the servicer excluding the reaction wheels

N = length(X0_serv)-6; % Number of Reaction Wheels
dt = Ttot(2)-Ttot(1); % Time step

X_serv = X0_serv; % Intializing the state vector for the servicer
X_deb = X0_deb; % Intializing the state vector for the target

SN = MRP2C(X_serv(1:3)); % DCM of the servicer from inertial frame to servicer body frame
DN = MRP2C(X_deb(1:3)); % DCM of the debris from inertial frame to debris body frame

params.attitude_mode = "First Attitude"; % Flag determining the attitude of the servicer
params.desat_flag = 0; % Flag determining if the servicer is in desat mode

params.sim.sig_RN = [0;0;0]; % Initialized reference attitude
params.sim.N_rvec_m = params.N_rvec_km*1000; % Initializing the seperation distance vector in m
params.sim.mode_code = 1; % 1: First, 2: Slew to Second, 3: Second, 4: Slew to First

results = storage(Ttot); % Allocating storage for the stored results

percent_check = 0; % Shows how much time is left in the simulations

params.sim.X_serv = X_serv; % store state of the servicer
params.sim.X_deb = X_deb; % store state of the target

% TODO: Not always true
% params.sim.H = [0,0,0];
for i_tt = 1:length(Ttot)
    
    t = Ttot(i_tt);
  
    
%%%% Desate Mode Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If the servicer is not in desat, check reaction wheel speeds
    RW_wheel_speeds = abs(X_serv(6+1:6+N,1));
    
    if params.desat_flag == 0
        for i_RW = 1:N
            if RW_wheel_speeds(i_RW) > params.wheel_speed_threshold
                params.desat_flag = 1;
                params.sat_wheel = i_RW;
            end
        end
    end
    
    % If the servicer is in desat, follow desat procedure
    if params.desat_flag == 1
        switch params.attitude_mode
            
            case "First Attitude"
%                 params.V = [10000,-10000];
                [data_anti,~,~] = find_anti_torque(SN'*params.sim.H,data);
                params.sim.sig_RN = C2MRP(data_anti.SN);
                
%                 params.sim.sig_RN = C2MRP(Euler3212C(deg2rad([-90,0,-20])));
                params.attitude_mode = "Slewing_To_Second_Attitude";
                
            case "Slewing_To_Second_Attitude"
%                 params.V = [0,0];
                RN = MRP2C(params.sim.sig_RN);
                SN = MRP2C(sig_SN);
                SR = SN*RN';
                sig_SR = C2MRP(SR);
                theta_SR = 4*atan(norm(sig_SR));
                % Once the spacecraft has reached the desired attitude
                if abs(theta_SR) < 0.035
                    params.attitude_mode = "Second Attitude";
                    params.wheel_speed_signs = sign(X_serv(6+1:6+N,1));
                end
                
            case "Second Attitude"
%                 params.V = [10000,-10000];
%                     if dot(SN'*params.sim.H/norm(SN'*params.sim.H),N_L_elect_serv/norm(N_L_elect_serv)) > 0 %as the target rotates, the torque changes which would take you out of this location without desat being finished
%                     if sum(sign(X_serv(6+1:6+N,1)) ~= sign(params.wheel_speed_signs)) > 0 %flips when one of the reaction wheels starts increasing in speed
%                    % if (normH-normHm1) > 0
%                         [data_anti,~,~] = find_anti_torque(SN'*params.sim.H,data);
%                         params.sim.sig_RN = C2MRP(data_anti.C2);
%                         %                 params.sim.sig_RN = [0.0911;-0.3205;-0.1814];
%                         params.attitude_mode = "Slewing_To_Third_Attitude";
%                     end
                    
                    if dot(SN'*params.sim.H/norm(SN'*params.sim.H),N_L_elect_serv/norm(N_L_elect_serv)) > 0 %as the target rotates, the torque changes which would take you out of this location without desat being finished
%                     if sign(X_serv(6+params.sat_wheel)) ~= sign(params.wheel_speed_signs(params.sat_wheel))
                        params.sim.sig_RN = [0;0;0];
                        params.attitude_mode = "Slewing_To_First_Attitude";
                    end
                
%             case "Slewing_To_Third_Attitude"
%                    %params.V = [0,0];
%                 RN = MRP2C(params.sim.sig_RN);
%                 SN = MRP2C(sig_SN);
%                 SR = SN*RN';
%                 sig_SR = C2MRP(SR);
%                 theta_SR = 4*atan(norm(sig_SR));
%                 % Once the spacecraft has reached the desired attitude
%                 if abs(theta_SR) < 0.0017
%                     params.attitude_mode = "Third Attitude";
%                     params.wheel_speed_signs = sign(X_serv(6+1:6+N,1));
%                 end
%                 
%             case "Third Attitude"
% %                 params.V = [10000,-10000];
%                     if dot(SN'*params.sim.H/norm(SN'*params.sim.H),N_L_elect_serv/norm(N_L_elect_serv)) > 0 %as the target rotates, the torque changes which would take you out of this location without desat being finished
%                         params.sim.sig_RN = [0;0;0];
%                         params.attitude_mode = "Slewing_To_First_Attitude";
%                     end
                
        
                
            case "Slewing_To_First_Attitude"
%                 params.V = [0,0];
                RN = MRP2C(params.sim.sig_RN);
                SN = MRP2C(sig_SN);
                SR = SN*RN';
                sig_SR = C2MRP(SR);
                theta_SR = 4*atan(norm(sig_SR));
                if abs(theta_SR) < 0.0017 
                    params.attitude_mode = "First Attitude";
                    params.desat_flag = 0;
%                     params.V = [10000,-10000];
                end
                
        end
        
    end

    switch params.attitude_mode
        case "First Attitude"
            params.sim.mode_code = 1;
        case "Slewing_To_Second_Attitude"
            params.sim.mode_code = 2;
        case "Second Attitude"
            params.sim.mode_code = 3;
        case "Slewing_To_First_Attitude"
            params.sim.mode_code = 4;
    end
    
%%%% Compute required control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    R_w_RN = [0,0,0]';
    R_w_RN_dot = [0,0,0]';
    
    sig_SN = X_serv(1:3);
    S_w_SN = X_serv(4:6);
    
    RN = MRP2C(params.sim.sig_RN);
    SN = MRP2C(sig_SN);
    SR = SN*RN';
    sig_SR = C2MRP(SR);

    S_w_RN = -SR*R_w_RN;
    S_w_SR = S_w_SN - S_w_RN;
    S_w_dot_RN = SR*R_w_RN_dot;
    
    Om_mat = X_serv(6+1:6+N,1);

    Gs = Gs0;
    hs = zeros(N,1);
    Hw = 0;
    for k = 1:N
    
        gsi = Gs(1:3,k);

        ws = dot(S_w_SN,gsi);

        hs(k) = Iws*(ws + Om_mat(k));
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    
    end
    
    
    
    [N_F_on_debris, N_F_on_serv, N_L_elect_debris, N_L_elect_serv, ~, overlapFlag] = ...
    multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.sim.N_rvec_m, params.V, DN', SN',...
    params.debris.D_COM, params.servicer.S_COM);
    
    % Convert each electrostatic torque into the corresponding body frame
    % used by each rigid-body EOM.
    S_L_serv = SN*N_L_elect_serv;
    D_L_deb = DN*N_L_elect_debris;  
    
    h_w = Gs*Iws*Om_mat;

    params.sim.Lr = (-K*sig_SR - P*(S_w_SR) - I_RW*((S_w_dot_RN) - tild(S_w_SN)*(S_w_RN)) + ...
    tild(S_w_SN)*(I_RW*S_w_SN + Gs*hs)-S_L_serv) + K_h*h_w;


    us = pinv(Gs)*-params.sim.Lr;
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


    params.sim.aterr = norm(sig_SR);
    params.sim.werr = norm(S_w_SR);

%%%% RK4 integrator for the servicing spacecraft %%%%%%%%%%%%%%%%%%%%%%%%%%

    Xdot_servicer = N_RW_EOM(X_serv, I_RW, Iws, Gs0, N,us,S_L_serv);

    k1 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k1/2, I_RW, Iws, Gs0, N,us,S_L_serv);
    
    k2 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k2/2, I_RW, Iws, Gs0, N,us,S_L_serv);

    k3 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k3, I_RW, Iws, Gs0, N,us,S_L_serv);
    k4 = Xdot_servicer*dt;

    X_serv = X_serv + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set if needed
    % TODO: Update so it switches the MRP based on the net forces
    % experienced during the rotation
    
    if norm(X_serv(1:3))>1
        X_serv(1:3) = -X_serv(1:3)/norm(X_serv(1:3))^2;
    end
    
        
    
    SN = MRP2C(X_serv(1:3));
    
%%%% RK4 integrator for the Target spacecraft %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xdot_deb = N_MRP_EOM(X_deb,params.debris.D_MI,D_L_deb);

    k1 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k1/2,params.debris.D_MI,D_L_deb);
    
    k2 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k2/2,params.debris.D_MI,D_L_deb);

    k3 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k3,params.debris.D_MI,D_L_deb);
    k4 = Xdot_deb*dt;

    X_deb = X_deb + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set
        if norm(X_deb(1:3))>1
           X_deb(1:3) = -X_deb(1:3)/norm(X_deb(1:3))^2;
        end

    DN = MRP2C(X_deb(1:3));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Store the updated state and momentum consistently after integration.
    params.sim.X_serv = X_serv;
    params.sim.X_deb = X_deb;

    S_w_SN = X_serv(4:6);
    w_DN = X_deb(4:6);
    Om_mat = X_serv(6+1:6+N,1);
    HB = I_RW*S_w_SN;
    Hw = zeros(3,1);
    for k = 1:N
        gsi = Gs0(1:3,k);
        ws = dot(S_w_SN,gsi);
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    end
    normHm1 = norm(params.sim.H);
    params.sim.H = HB + Hw;
    normH = norm(params.sim.H);
    H_deb = params.debris.D_MI*w_DN;
    
    params.sim.dot_prod = dot(N_L_elect_serv/norm(N_L_elect_serv),(SN'*params.sim.H)/norm(SN'*params.sim.H));

    params.sim.H_deb = H_deb;
    
    params.sim.us = us;
    params.sim.L_e_serv = S_L_serv;
    params.sim.L_e_deb = D_L_deb;
 
    results = update_storage(results,i_tt,params);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Real Time Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_left_percentage = t/Ttot(end);
    if time_left_percentage*100 > percent_check || t == Ttot(end)

        disp("Time  "+ t +"s:   "+time_left_percentage*100 + "% complete")
        percent_check = percent_check + 1;
    
        
        N_H = SN'*params.sim.H;
        
        N_servicer_COM = SN'*params.servicer.S_COM + params.N_rvec_km*1000;
        N_debris_COM = DN'*params.debris.D_COM;
        
        fig = figure(100);
        clf(fig)
        set(fig,'Position',[50 50 2000 900])
        
        subplot(2,6,[1,2,7,8])
        hold on
        set(gca,'FontName','times')
        makeSphsPicture_2craft(params.debris.D_spheres, params.servicer.S_spheres,...
            [0,0,0],params.N_rvec_km*1000, [params.debris.voltage, params.servicer.voltage], DN', SN')

        quiver3(N_servicer_COM(1),N_servicer_COM(2),N_servicer_COM(3),...
            10000*N_L_elect_serv(1),10000*N_L_elect_serv(2),10000*N_L_elect_serv(3),'Linewidth',2)
        
        quiver3(N_servicer_COM(1),N_servicer_COM(2),N_servicer_COM(3),...
            10*N_H(1),10*N_H(2),10*N_H(3),'Linewidth',2)
        
        
%         quiver3(S_servicer_COM(1)+params.N_rvec_km(1)*1000,S_servicer_COM(2)+params.N_rvec_km(2)*1000,S_servicer_COM(3)+params.N_rvec_km(3)*1000,...
%             10000*N_F_on_serv(1),10000*N_F_on_serv(2),10000*N_F_on_serv(3),'Linewidth',2)
%         
%         quiver3(D_debris_COM(1),D_debris_COM(2),D_debris_COM(3),...
%             10000*N_L_elect_debris(1),10000*N_L_elect_debris(2),10000*N_L_elect_debris(3),'Linewidth',2)
%         
%         quiver3(D_debris_COM(1),D_debris_COM(2),D_debris_COM(3),...
%             10000*N_F_on_debris(1),10000*N_F_on_debris(2),10000*N_F_on_debris(3),'Linewidth',2)
%         
        scatter3(N_debris_COM(1),N_debris_COM(2),N_debris_COM(3),20,'k','filled')
        scatter3(N_servicer_COM(1),N_servicer_COM(2),N_servicer_COM(3),20,'k','filled')

        
        axis equal
        xlim([-20,50])
        ylim([-20,20])
        zlim([-20,20])
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Z [m]')
        title('Spacecraft Geometry')
        
        view(3)
        hold off

        subplot(2,6,[3,4])
        H_hist = vecnorm(results.Htot_servicer(:,1:i_tt),2,1);
        plot(results.Ttot(1:i_tt)/3600,H_hist,'LineWidth',2)
        box on
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('H (Nms)')
        title('Angular Momentum History')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
        
        subplot(2,6,[5,6])
        hold on
        rw_hist = results.Xtot_servicer(7:9,1:i_tt);
        plot(results.Ttot(1:i_tt)/3600,rw_hist(1,:),'LineWidth',2)
        plot(results.Ttot(1:i_tt)/3600,rw_hist(2,:),'LineWidth',2)
        plot(results.Ttot(1:i_tt)/3600,rw_hist(3,:),'LineWidth',2)
        yline(params.wheel_speed_threshold,'--k')
        yline(-params.wheel_speed_threshold,'--k')
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('Wheel Speed (rad/s)')
        title('Servicer Reaction Wheel Speeds')
        legend('$\Omega_1$','$\Omega_2$','$\Omega_3$','Location','best')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
%         ylim([-1.1*params.wheel_speed_threshold,params.wheel_speed_threshold*1.1])
        
        subplot(2,6,[9,10])
        hold on
        dot_prod_hist = results.dot_prods(:,1:i_tt);
        plot(results.Ttot(1:i_tt)/3600,dot_prod_hist,'LineWidth',2)
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('dot($^BH_{serv}$,$^BL_{serv})$')
        title('H $L_e$ angle')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
        ylim([-1,1])
        
        subplot(2,6,[11,12])
        hold on
        plot(results.Ttot(1:i_tt)/3600,results.debris_ang_vel(1,1:i_tt)*180/pi,'LineWidth',2)
        plot(results.Ttot(1:i_tt)/3600,results.debris_ang_vel(2,1:i_tt)*180/pi,'LineWidth',2)
        plot(results.Ttot(1:i_tt)/3600,results.debris_ang_vel(3,1:i_tt)*180/pi,'LineWidth',2)
        plot(results.Ttot(1:i_tt)/3600,results.debris_ang_speed(1:i_tt)*180/pi,'k','LineWidth',2)
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('$\omega_{DN}$ (deg/s)')
        title('Debris Angular Velocity')
        legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','best')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
        
    end


end
end
