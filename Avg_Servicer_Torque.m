% Finding the Average Torque the servicer experiences at a give attitude
%%% average based on all possible servicer attitudes
%%% this is for finding servicer attitudes that generate consistent torques
%%% regardless of target attitude


function[data, Ls_norms, overlap_flag] = Avg_Servicer_Torque(params,n,SN)

[all_torques, Ls_norms, ~, overlap_flag] = All_Torques_Target(params,n,SN);

N_F_on_debris_all_target_atts = zeros(n^3,3);
N_F_on_serv_all_target_atts = zeros(n^3,3);
N_L_elect_debris_all_target_atts = zeros(n^3,3);
N_L_elect_serv_all_target_atts = zeros(n^3,3);

for i_D_EAs = 1:n^3
    N_F_on_debris_all_target_atts(i_D_EAs,:) = all_torques{i_D_EAs}.N_F_on_debris;
    N_F_on_serv_all_target_atts(i_D_EAs,:) = all_torques{i_D_EAs}.N_F_on_serv;
    N_L_elect_debris_all_target_atts(i_D_EAs,:) = all_torques{i_D_EAs}.N_L_elect_debris;
    N_L_elect_serv_all_target_atts(i_D_EAs,:) = all_torques{i_D_EAs}.N_L_elect_serv;
end
data.avg_N_F_on_debris = mean(N_F_on_debris_all_target_atts);
data.avg_N_F_on_serv = mean(N_F_on_serv_all_target_atts);
data.avg_N_L_elect_debris = mean(N_L_elect_debris_all_target_atts);
data.avg_N_L_elect_serv = mean(N_L_elect_serv_all_target_atts);
data.overlapFlag = overlap_flag;
data.all_Ls_target_rotating = all_torques;

data.Ls_target_rotation_norm = mean(Ls_norms);

end









