function results = N_RW_sim(X0_serv,X0_deb, Iws, Gs0, K, P, Ttot,params,data)

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

for i_tt = 1:length(Ttot)
  
    t = Ttot(i_tt); % Determine current time
    if t == 50*3600
        ddd = 0;
    end
    
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
%                 [data_anti,~,~] = find_anti_torque(B_L2,data);
%                 params.sim.sig_RN = C2MRP(data_anti.C2);
                params.sim.sig_RN = [0,0,1]';
                params.attitude_mode = "Slewing_To_Second_Attitude";
                
            case "Slewing_To_Second_Attitude"
%                 params.V = [0,0];
                RN = MRP2C(params.sim.sig_RN);
                SN = MRP2C(sig_SN);
                SR = SN*RN';
                sig_SR = C2MRP(SR);
                theta_SR = 4*atan(norm(sig_SR));
                % Once the spacecraft has reached the desired attitude
                if abs(theta_SR) < 0.0017
                    params.attitude_mode = "Second Attitude";
                    params.wheel_speed_signs = sign(X_serv(6+1:6+N,1));
                end
                
            case "Second Attitude"
%                 params.V = [10000,-10000];
%                 if dot(X_serv(7:9)/norm(X_serv(7:9)),L_serv/norm(L_serv)) > 0 as the target rotates, the torque changes which would take you out of this location without desat being finished
                 if sum(sign(X_serv(6+1:6+N,1)) ~= sign(params.wheel_speed_signs)) > 0
                    params.sim.sig_RN = [0;0;0];
                    params.attitude_mode = "Slewing_To_First_Attitude";
                end
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
                    params.V = [10000,-10000];
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
    
    [D_F_on_debris, S_F_on_serv, D_L_elect_debris, S_L_elect_serv, ~, overlapFlag] = ...
    multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.sim.N_rvec_m, params.V, DN, SN,...
    params.debris.D_COM, params.servicer.S_COM);
    
    % Convert each electrostatic torque into the corresponding body frame
    % used by each rigid-body EOM.
    L_serv = S_L_elect_serv;
    L_deb = D_L_elect_debris;        

    params.sim.Lr = -K*sig_SR - P*(S_w_SR) + I_RW*((S_w_dot_RN) - tild(S_w_SN)*(S_w_RN)) + ...
    tild(S_w_SN)*(I_RW*S_w_SN + Gs*hs)-L_serv;

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

    Xdot_servicer = N_RW_EOM(X_serv, I_RW, Iws, Gs0, N,us,L_serv);

    k1 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k1/2, I_RW, Iws, Gs0, N,us,L_serv);
    
    k2 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k2/2, I_RW, Iws, Gs0, N,us,L_serv);

    k3 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_serv+k3, I_RW, Iws, Gs0, N,us,L_serv);
    k4 = Xdot_servicer*dt;

    X_serv = X_serv + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set if needed
        if norm(X_serv(1:3))>1
           X_serv(1:3) = -X_serv(1:3)/norm(X_serv(1:3))^2;
        end
    
    SN = MRP2C(X_serv(1:3));
    
%%%% RK4 integrator for the Target spacecraft %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xdot_deb = N_MRP_EOM(X_deb,params.debris.D_MI,L_deb);

    k1 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k1/2,params.debris.D_MI,L_deb);
    
    k2 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k2/2,params.debris.D_MI,L_deb);

    k3 = Xdot_deb*dt;
    Xdot_deb = N_MRP_EOM(X_deb+k3,params.debris.D_MI,L_deb);
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
    
    params.sim.H = HB + Hw;
    
    H_deb = params.debris.D_MI*w_DN;

    params.sim.H_deb = H_deb;
    
    params.sim.us = us;
    params.sim.L_e_serv = L_serv;
    params.sim.L_e_deb = L_deb;
    
 
    results = update_storage(results,i_tt,params);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Real Time Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%
    time_left_percentage = t/Ttot(end);
    if time_left_percentage*100 > percent_check

        disp("Time  "+ t +"s:   "+time_left_percentage*100 + "% complete")
        percent_check = percent_check +2;
      
        plotting_servicer = params.servicer.N_spheres;
        for i = 1:length(params.servicer.N_spheres)
            sph_loc = params.servicer.N_spheres(1:3,i);
            new_sph_loc_serv = SN*sph_loc;
            plotting_servicer(1:3,i) = new_sph_loc_serv;
        end
        
        plotting_debris = params.debris.N_spheres;
        for i = 1:length(params.debris.N_spheres)
            sph_loc = params.debris.N_spheres(1:3,i);
            new_sph_loc_deb = DN*sph_loc;
            plotting_debris(1:3,i) = new_sph_loc_deb;
        end
        
        N_servicer_COM = SN'*params.servicer.S_COM;
        N_debris_COM = DN'*params.debris.D_COM;
        
        fig = figure(100);
        clf(fig)
        set(fig,'Position',[50 50 1000 500])
        
        subplot(2,4,[1,2,5,6])
        hold on
        set(gca,'FontName','times')
        makeSphsPicture_2craft(plotting_debris, plotting_servicer,...
            N_debris_COM,N_servicer_COM+params.N_rvec_km*1000, [params.debris.voltage, params.servicer.voltage])

        quiver3(N_servicer_COM(1)+params.N_rvec_km(1)*1000,N_servicer_COM(2)+params.N_rvec_km(2)*1000,N_servicer_COM(3)+params.N_rvec_km(3)*1000,...
            10000*L_serv(1),10000*L_serv(2),10000*L_serv(3),'Linewidth',2)
        
        quiver3(N_servicer_COM(1)+params.N_rvec_km(1)*1000,N_servicer_COM(2)+params.N_rvec_km(2)*1000,N_servicer_COM(3)+params.N_rvec_km(3)*1000,...
            10000*S_F_on_serv(1),10000*S_F_on_serv(2),10000*S_F_on_serv(3),'Linewidth',2)
        
        quiver3(N_debris_COM(1),N_debris_COM(2),N_debris_COM(3),...
            10000*D_L_elect_debris(1),10000*D_L_elect_debris(2),10000*D_L_elect_debris(3),'Linewidth',2)
        
        quiver3(N_debris_COM(1),N_debris_COM(2),N_debris_COM(3),...
            10000*D_F_on_debris(1),10000*D_F_on_debris(2),10000*D_F_on_debris(3),'Linewidth',2)
        
        axis equal
        xlim([-3,50])
        ylim([-17,17])
        zlim([-12,12])
%         c=colorbar;
%         c.Label.String = 'Surface Charge Density (nC/m^2)';
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Z [m]')
        title('Spacecraft Geometry')
        
        view(3)
        hold off

        subplot(2,4,[3,4])
        H_hist = vecnorm(results.Htot_servicer(:,1:i_tt),2,1);
        plot(results.Ttot(1:i_tt)/3600,H_hist,'LineWidth',2)
        box on
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('H (Nms)')
        title('Angular Momentum History')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
        
        subplot(2,4,[7,8])
        debris_w_hist = vecnorm(results.debris_ang_vel(:,1:i_tt),2,1);
        plot(results.Ttot(1:i_tt)/3600,debris_w_hist,'LineWidth',2)
        set(gca,'FontName','times')
        xlabel('Time (hours)')
        ylabel('$\omega_{DN}$ (rad/s)')
        title('Debris Angular Velocity')
        xlim([results.Ttot(1), results.Ttot(end)]/3600)
        
    end


end
end
