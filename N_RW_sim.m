function results = N_RW_sim(X0_servicer,X0_target, Iws, Gs0, K, P, Ttot,params,data)

I_RW = params.servicer.S_MI; % The moment of inertia of the servicer excluding the reaction wheels

N = length(X0_servicer)-6; % Number of Reaction Wheels
dt = Ttot(2)-Ttot(1); % Time step

X_servicer = X0_servicer; % Intializing the state vector for the servicer
X_target = X0_target; % Intializing the state vector for the target

SN = MRP2C(X_servicer(1:3)); % DCM of the servicer from inertial frame to servicer body frame
DN = MRP2C(X_target(1:3)); % DCM of the debris from inertial frame to debris body frame

params.attitude_mode = "First Attitude"; % Flag determining the attitude of the servicer
params.desat_flag = 0; % Flag determining if the servicer is in desat mode

params.sim.sig_RN = [0;0;0]; % Initialized reference attitude
params.sim.r_m = params.N_rvec_km*1000; % Initializing the seperation distance vector in m
params.sim.mode_code = 1; % 1: First, 2: Slew to Second, 3: Second, 4: Slew to First

results = storage(Ttot); % Allocating storage for the stored results

percent_check = 0; % Shows how much time is left in the simulations

params.sim.X_servicer = X_servicer; % store state of the servicer
params.sim.X_target = X_target; % store state of the target

for i_tt = 1:length(Ttot)
  
    t = Ttot(i_tt); % Determine current time
    if t == 7895
        ddd = 0;
    end
    
%%%% Desate Mode Check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If the servicer is not in desat, check reaction wheel speeds
    RW_wheel_speeds = abs(X_servicer(6+1:6+N,1));
    
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
                params.V = [10000,-10000];
%                 [data_anti,~,~] = find_anti_torque(B_L2,data);
%                 params.sim.sig_RN = C2MRP(data_anti.C2);
                params.sim.sig_RN = [0,0,1]';
                params.attitude_mode = "Slewing_To_Second_Attitude";
                
            case "Slewing_To_Second_Attitude"
                params.V = [0,0];
                RN = MRP2C(params.sim.sig_RN);
                SN = MRP2C(sig_SN);
                SR = SN*RN';
                sig_SR = C2MRP(SR);
                theta_SR = 4*atan(norm(sig_SR));
                % Once the spacecraft has reached the desired attitude
                if abs(theta_SR) < 0.0017
                    params.attitude_mode = "Second Attitude";
                end
                
            case "Second Attitude"
                params.V = [10000,-10000];
                if dot(X_servicer(7:9)/norm(X_servicer(7:9)),B_L2/norm(B_L2)) > 0
                    params.sim.sig_RN = [0;0;0];
                    params.attitude_mode = "Slewing_To_First_Attitude";
                end
            case "Slewing_To_First_Attitude"
                params.V = [0,0];
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
    
    sig_SN = X_servicer(1:3);
    w_SN = X_servicer(4:6);
    
    RN = MRP2C(params.sim.sig_RN);
    SN = MRP2C(sig_SN);
    SR = SN*RN';
    sig_SR = C2MRP(SR);

    S_w_RN = -SR*R_w_RN;
    S_w_SR = w_SN - S_w_RN;
    w_dot_RN = SR*R_w_RN_dot;
    
    Om_mat = X_servicer(6+1:6+N,1);

    Gs = Gs0;
    hs = [];
    Hw = 0;
    for k = 1:N
    
        gsi = Gs(1:3,k);

        ws = dot(w_SN,gsi);

        hsi = Iws*(ws + Om_mat(k));
        hs = [hs;hsi];
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    
    end
    
    [D_F_on_debris, S_F_on_serv, D_L_elect_debris, S_L_elect_serv, ~, overlapFlag] = ...
    multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.sim.r_m, params.V, DN, SN,...
    params.debris.D_COM, params.servicer.S_COM);
    
    % Convert each electrostatic torque into the corresponding body frame
    % used by each rigid-body EOM.
    L_serv = S_L_elect_serv;
    L_debris = D_L_elect_debris;
    B_L2 = L_serv;
        

    params.sim.Lr = -K*sig_SR - P*(S_w_SR) + I_RW*((w_dot_RN) - tild(w_SN)*(S_w_RN)) + ...
    tild(w_SN)*(I_RW*w_SN + Gs*hs)-L_serv;

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

    Xdot_servicer = N_RW_EOM(X_servicer, I_RW, Iws, Gs0, N,us,L_serv);

    k1 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_servicer+k1/2, I_RW, Iws, Gs0, N,us,L_serv);
    
    k2 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_servicer+k2/2, I_RW, Iws, Gs0, N,us,L_serv);

    k3 = Xdot_servicer*dt;
    Xdot_servicer = N_RW_EOM(X_servicer+k3, I_RW, Iws, Gs0, N,us,L_serv);
    k4 = Xdot_servicer*dt;

    X_servicer = X_servicer + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set if needed
        if norm(X_servicer(1:3))>1
           X_servicer(1:3) = -X_servicer(1:3)/norm(X_servicer(1:3))^2;
        end
    
    SN = MRP2C(X_servicer(1:3));
    
%%%% RK4 integrator for the Target spacecraft %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Xdot_target = N_MRP_EOM(X_target,params.debris.D_MI,L_debris);

    k1 = Xdot_target*dt;
    Xdot_target = N_MRP_EOM(X_target+k1/2,params.debris.D_MI,L_debris);
    
    k2 = Xdot_target*dt;
    Xdot_target = N_MRP_EOM(X_target+k2/2,params.debris.D_MI,L_debris);

    k3 = Xdot_target*dt;
    Xdot_target = N_MRP_EOM(X_target+k3,params.debris.D_MI,L_debris);
    k4 = Xdot_target*dt;

    X_target = X_target + 1/6*(k1+2*k2+2*k3+k4);
    
    %Switch to shadow set
        if norm(X_target(1:3))>1
           X_target(1:3) = -X_target(1:3)/norm(X_target(1:3))^2;
        end

    DN = MRP2C(X_target(1:3));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Store the updated state and momentum consistently after integration.
    params.sim.X_servicer = X_servicer;
    params.sim.X_target = X_target;

    w_SN = X_servicer(4:6);
    Om_mat = X_servicer(6+1:6+N,1);
    HB = I_RW*w_SN;
    Hw = zeros(3,1);
    for k = 1:N
        gsi = Gs0(1:3,k);
        ws = dot(w_SN,gsi);
        Hw = Hw + Iws*(ws + Om_mat(k))*gsi;
    end
    
    params.sim.H = HB + Hw;
    
    params.sim.us = us;
    params.sim.L_e = L_serv;
    
 
    results = update_storage(results,i_tt,params);
    
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
        c=colorbar;
        c.Label.String = 'Surface Charge Density (nC/m^2)';
        xlabel('X [m]')
        ylabel('Y [m]')
        zlabel('Z [m]')
        
        view(3)
        hold off
        
        
    end


end
end
