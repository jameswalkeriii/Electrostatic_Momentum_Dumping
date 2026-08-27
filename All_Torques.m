% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function [data,Ls,EAs] = All_Torques(params,n)
n_p = round(n/2);
range_yaw = linspace(-pi,pi,n);
range_pitch = linspace(-pi/2,pi/2,n_p);
range_roll = linspace(-pi,pi,n);
i = 1;
Ls = zeros(n^2*n_p,3);
data = cell(1,n^2*n_p);
EAs = zeros(n^2*n_p,3);

ND = eye(3);
flag = 0;
for y_s = 1:n
    yaw_servicer = range_yaw(y_s);
    for p_s = 1:n_p
        pitch_servicer = range_pitch(p_s);
        for r_s = 1:n
            roll_servicer = range_roll(r_s);
            EA_servicer = [yaw_servicer,pitch_servicer,roll_servicer];
            SN = Euler3212C(EA_servicer);
            NS = SN';
            
            [data_target_rotations, Ls_norms, ~,all_torques_target_rotation] = Avg_Servicer_Torque(params,n,SN);
            
            avg_Ls_target_rotation_norm = mean(Ls_norms);
            avg_Ls_target_rotation_norm = avg_Ls_target_rotation_norm/norm(avg_Ls_target_rotation_norm);
            
            
            
            rot_angs = zeros(n^2*n_p,1);
            for i_targets = 1:length(Ls_norms)
                L_tr = all_torques_target_rotation{i_targets}.N_L_elect_serv;
                L_tr_norm = L_tr/norm(L_tr);
                
                %                 data_target_rotations.all_Ls_target_rotating{i_targets}.angle_from_mean = rad2deg(acos(dot(avg_Ls_target_rotation_norm,L_tr_norm)));
                angle_from_mean = rad2deg(acos(dot(avg_Ls_target_rotation_norm,L_tr_norm)));
                rot_angs(i_targets) = angle_from_mean;
                if angle_from_mean > 80
                    flag = 1;
                end
                    
            end
            
            average_dist_angle = mean(rot_angs);
            std_dist_angle = std(rot_angs);
            
            
            if flag == 0
                % 0 = narrow distribution
                data{i} = data_target_rotations;
                data{i}.dist = 0;
                data{i}.EA = EA_servicer;
                data{i}.NS = NS;
                data{i}.SN = SN;
                data{i}.MRP = C2MRP(SN);
                Ls(i,:) = avg_Ls_target_rotation_norm;
                EAs(i,:) = EA_servicer;
                data{i}.average_dist_angle = average_dist_angle;
                data{i}.std_dist_angle = std_dist_angle;
                if average_dist_angle > 50
                    data{i}.dist = 2;
                end
            else
                % 1 = narrow distribution
                data{i} = data_target_rotations;
                data{i}.dist = 1;
                data{i}.EA = EA_servicer;
                data{i}.NS = NS;
                data{i}.SN = SN;
                data{i}.MRP = C2MRP(SN);
                Ls(i,:) = avg_Ls_target_rotation_norm;
                EAs(i,:) = EA_servicer;
                data{i}.average_dist_angle = average_dist_angle;
                data{i}.std_dist_angle = std_dist_angle;
                
            end
            flag = 0;
            
            i = i+1;
        end
    end
end

data = data(~cellfun('isempty', data));

end
