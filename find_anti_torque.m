function [data_anti,i_anti,tot] = find_anti_torque(H,data)
%Must imput a non normalized H
% MAKE SURE ALL IN THE SAME FRAME
anti_parr_check_min = -1;
data_anti = 0;
i_anti = 0;
tot = zeros(length(data),1);
dot_product_norm = zeros(length(data),1);
dot_product_magnitude = zeros(length(data),1);

for i = 1:length(data)
    % If the spacecraft are overlapping, then set the value as NaN
    if data{i}.overlapFlag == 0
        % If not, add the sum measure to tot
        E_torque_Nframe = data{i}.N_L_elect_serv;
        Ang_mom_Nframe = H;
        
        E_torque_Nframe_dir = E_torque_Nframe/norm(E_torque_Nframe);
        Ang_mom_Nframe_dir = Ang_mom_Nframe/norm(Ang_mom_Nframe);
        
        dot_product_norm(i) = -dot(E_torque_Nframe_dir,Ang_mom_Nframe_dir);
        
        dot_product_magnitude(i) = norm(dot(E_torque_Nframe,Ang_mom_Nframe));
    else
        dot_product_norm(i) = NaN;
        dot_product_magnitude(i) = NaN;
    end
end

dot_product_magnitude_percent = dot_product_magnitude/max(dot_product_magnitude);

for i = 1:length(data)
    
    weighted_value = 10*dot_product_norm(i)+ dot_product_magnitude_percent(i);
    
    tot(i,1:4) = [dot_product_norm(i),dot_product_magnitude(i),dot_product_magnitude_percent(i),weighted_value];
    if weighted_value > anti_parr_check_min
        clear data_anti
        data_anti = data{i};
        i_anti = i;
        anti_parr_check_min = weighted_value;
    end
    
end
