% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function [data,Ls,EAs] = All_Torques(params,n)
range_yaw = linspace(-pi,pi,n);
range_pitch = linspace(-pi/2,pi/2,n);
range_roll = linspace(-pi,pi,n);
i = 1;
Ls = zeros(n^3,3);
data = cell(1,n^3);
EAs = zeros(n^3,3);

ND = eye(3);

for y_s = 1:n
    yaw_servicer = range_yaw(y_s);
    for p_s = 1:n
        pitch_servicer = range_pitch(p_s);
        for r_s = 1:n
            roll_servicer = range_roll(r_s);
            EA_servicer = [yaw_servicer,pitch_servicer,roll_servicer];
            SN = Euler3212C(EA_servicer);
            NS = SN';
            
            [data_target_rotations, Ls_norms] = Avg_Servicer_Torque(params,n,SN);

            avg_Ls_target_rotation_norm = mean(Ls_norms);
            
            data{i} = data_target_rotations;
            data{i}.EA = EA_servicer;
            data{i}.NS = NS;
            data{i}.SN = SN;
            data{i}.MRP = C2MRP(SN);
            Ls(i,:) = avg_Ls_target_rotation_norm;
%             sigs(i,:) = C2MRP(Euler3212C(EA));
            EAs(i,:) = EA_servicer;
            
            for i_targets = 1:length(Ls_norms)
                L_tr = data_target_rotations.all_Ls_target_rotating{i_targets}.N_L_elect_serv;
                L_tr_norm = L_tr/norm(L_tr);
                
                data_target_rotations.all_Ls_target_rotating{i_targets}.angle_from_mean = rad2deg(acos(dot(avg_Ls_target_rotation_norm,L_tr_norm)));
                rot_angs(i_targets) = data_target_rotations.all_Ls_target_rotating{i_targets}.angle_from_mean;
            end
            
            data{i}.average_dist_angle = mean(rot_angs);
            
            if data{i}.average_dist_angle > 40
                data{i}.dist = 1;
                % 1 = wide distribution
            else
                data{i}.dist = 0;
                % 0 = narrow distribution
            end
            
            i = i+1;
        end
    end
end

end
