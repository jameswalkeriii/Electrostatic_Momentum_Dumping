function [Xtot, Htot, Ttot, ttot,aterrtot,werrtot,ustot,Lrtot,...
    L_e_tot,params,normBN,reftot,rtot,flags] = N_RW_sim(X0, I_RW, Iws, Gs0,...
    N, K, P, tn, dt,params,data)

% Allocating parameters to respective variables
I_RW = params.servicer.B_MI;


X = X0; 
ttot = 0; 
Xtot = X0; Htot = []; Ttot = []; aterrtot = []; werrtot = []; ustot = [];
Lrtot = []; L_e_tot = []; normBN = []; reftot = []; rtot = []; flags = [];

r = params.r_km*1000; V = params.V; C1 = eye(3); C2 = MRP2C(X(1:3));

params.attitude_mode = "First Attitude";
params.desat_mode = "IDLE";
params.desat_flag = 0;

r0 = r;
t = 0; 
refflag = 0; wheel_speed_flag = 0;
sig_RN = [0;0;0]; sig_anti = [ -1;0;0];
params.sig_RN = sig_RN;

% sig_anti = [-0.0186;-0.9605;-0.0297];
Xprv = X0;
percent_check = 0;

[ ~, ~, ~, L20, ~, overlapFlag] = ...
    multisphereFT( params.debris.spheres, params.servicer.spheres, r, V, C1, C2,...
    params.debris.COM, params.servicer.COM);
L = L20;
% L2 =  [ -0.000054008811829;0.003719000304351;-0.002708315879101];
while t < tn
    
    if params.desat_flag == 0
        for i_RW = 1:N
            if abs(X(6+i_RW)) > params.wheel_speed_threshold
                params.desat_flag = 1;
                params.sat_wheel = i_RW;
            end
        end
    end
    
    if params.desat_flag == 1
        switch params.attitude_mode
            
            case "First Attitude"
                [data_anti,~,~] = find_anti_torque(B_L2,data);
                params.sig_RN = C2MRP(data_anti.C2);
                params.attitude_mode = "Slewing_To_Second_Attitude";
                
            case "Slewing_To_Second_Attitude"
                if abs(norm(X(1:3)) - norm(params.sig_RN)) < .001 % TODO: Look at angle between not just norm
                    params.attitude_mode = "Second Attitude";
                end
                
            case "Second Attitude"
                if dot(X(7:9)/norm(X(7:9)),B_L2/norm(B_L2)) > 0
                    params.sig_RN = [0;0;0];
                    params.attitude_mode = "Slewing_To_First_Attitude";
                end
            case "Slewing_To_First_Attitude"
                if abs(norm(X(1:3)) - norm(params.sig_RN)) < .001 % TODO: Look at angle between not just norm
                    params.attitude_mode = "First Attitude";
                    params.desat_flag = 0;
                end
                
        end
        
    end

    

    sig_RN = params.sig_RN;
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
    [ F1, F2, T_L1, T_L2, ~, overlapFlag] = ...
    multisphereFT( params.debris.spheres, params.servicer.spheres, r, V, C1, C2,...
    params.debris.COM, params.servicer.COM);
    
    B_L2 = C2*T_L2;

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
    
    time_left_percentage = t/tn;
    if time_left_percentage*100 > percent_check
        disp("Time  "+ t +"s:   "+time_left_percentage*100 + "% complete")
        percent_check = percent_check +1;
      
        plotting_servicer = params.servicer.spheres;
        for i = 1:length(params.servicer.spheres)
            sph_loc = params.servicer.spheres(1:3,i);
            new_sph_loc = C2*sph_loc;
            plotting_servicer(1:3,i) = new_sph_loc;
        end
        
        fig = figure(100);
        clf(fig)
        hold on
        set(gca,'FontName','times')
        makeSphsPicture_2craft(plotting_servicer, params.debris.spheres, params.r_km*1000,...
            [0 0 0], [params.servicer.voltage, params.debris.voltage])
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

