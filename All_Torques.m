% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function [data,Ls,sigs,EAs] = All_Torques(params,n)
range_yaw = linspace(-pi,pi,n);
range_pitch = linspace(-pi/2,pi/2,n);
range_roll = linspace(-pi,pi,n);
i = 1;
Ls = zeros(n^3,3);
data = cell(1,n^3);
EAs = zeros(n^3,3);

for y = 1:n
    yaw = range_yaw(y);
    for p = 1:n
        pitch = range_pitch(p);
        for r = 1:n
            roll = range_roll(r);
            EA = [yaw,pitch,roll];
            SN = Euler3212C(EA);
            [N_F_on_debris, N_F_on_serv, N_L_elect_debris, N_L_elect_serv, qs, overlapFlag] = ...
                multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.N_rvec_km*1000,...
                params.V, eye(3), SN, params.debris.D_COM, SN'*params.servicer.S_COM);
            data{i}.EA = EA;
            data{i}.C2 = SN;
            data{i}.MRP = C2MRP(SN');
            data{i}.N_F_on_debris = N_F_on_debris;
            data{i}.N_F_on_serv = N_F_on_serv;
            data{i}.N_L_elect_debris = N_L_elect_debris;
            data{i}.N_L_elect_serv = N_L_elect_serv;
            data{i}.overlapFlag = overlapFlag;
            Ls(i,:) = N_L_elect_serv/norm(N_L_elect_serv);
            sigs(i,:) = C2MRP(Euler3212C(EA)');
            EAs(i,:) = EA;
            i = i+1;
        end
    end
end

end
