function [data_anti,i_anti,tot] = find_anti_torque(H,data)
%Must imput a non normalized H
% MAKE SURE ALL IN THE SAME FRAME
anti_parr_check_min = 1;
data_anti = 0;
i_anti = 0;
tot = zeros(length(data),4);

for i = 1:length(data)
    % If the spacecraft are overlapping, then set the value as NaN
    if data{i}.overlapFlag == 0
        % If not, add the sum measure to tot
        E_torque_Bframe = data{i}.B2_L2;
        Ang_mom_Bframe = H;
        
        E_torque_Bframe_dir = E_torque_Bframe/norm(E_torque_Bframe);
        Ang_mom_Bframe_dir = Ang_mom_Bframe/norm(Ang_mom_Bframe);
        
        tot(i) = dot(E_torque_Bframe_dir,Ang_mom_Bframe_dir);
        
        
        if tot(i) < anti_parr_check_min
            clear data_anti
            data_anti = data{i};
            i_anti = i;
            anti_parr_check_min = tot(i);
        end
    else
        tot(i,1:4) = NaN;
    end
    
    if anti_parr_check_min > 0
        disp('No torques acting opposite the angular momentum')
    end
end