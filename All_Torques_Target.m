% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function [data,Ls,EAs,overlap_flag] = All_Torques_Target(params,n,SN)
n_p = round(n/2);
range_yaw = linspace(-pi,pi,n);
range_pitch = linspace(-pi/2,pi/2,n_p);
range_roll = linspace(-pi,pi,n);
i = 1;
Ls = zeros(n^2*n_p,3);
data = cell(1,n^2*n_p);
EAs = zeros(n^2*n_p,3);

overlap = 0;

NS = SN';

for y = 1:n
    yaw = range_yaw(y);
    for p = 1:n_p
        pitch = range_pitch(p);
        for r = 1:n
            roll = range_roll(r);
            EA = [yaw,pitch,roll];
%             if EA(1) == 0
%                 if EA(2) == 0 
%                     if EA(3) > pi-.001
%                         ddd=0;
%                     elseif EA(3) >0
%                         ddd = 0;
%                     end
%                 end
%             end
            DN = Euler3212C(EA);
            ND = DN';
            [N_F_on_debris, N_F_on_serv, N_L_elect_debris, N_L_elect_serv, qs, overlapFlag] = ...
                multisphereFT( params.debris.D_spheres, params.servicer.S_spheres, params.N_rvec_km*1000,...
                params.V, ND, NS, params.debris.D_COM, params.servicer.S_COM);
            data{i}.EA = EA;
            data{i}.ND = ND;
            data{i}.DN = DN;
            data{i}.MRP = C2MRP(DN);
            data{i}.N_F_on_debris = N_F_on_debris;
            data{i}.N_F_on_serv = N_F_on_serv;
            data{i}.N_L_elect_debris = N_L_elect_debris;
            data{i}.N_L_elect_serv = N_L_elect_serv;
            data{i}.overlapFlag = overlapFlag;
            Ls(i,:) = N_L_elect_serv/norm(N_L_elect_serv);
%             sigs(i,:) = C2MRP(Euler3212C(EA));
            EAs(i,:) = EA;
            i = i+1;
            overlap = overlap+overlapFlag;

%             PlotInitialPosition(params, SN, DN, N_L_elect_serv);
            
            
        end
    end
end

overlap_flag = sum(overlap);

end
