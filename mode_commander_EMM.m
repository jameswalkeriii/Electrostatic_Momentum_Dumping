function params = mode_commander_EMM(params)






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
%                 [data_anti,~,~] = find_anti_torque(L,data);
                sig_RN = sig_anti; 
                refflag = 1;
                sat_wheel = (i+6);
                i = 4;
%                 X(1:3) = sig_anti;
%                 sig_RN = sig_anti;
%                 X(1:3) = [0,0,0]';
%                 sig_RN = [0,0,0]';
                r(1) = 10000;
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



end