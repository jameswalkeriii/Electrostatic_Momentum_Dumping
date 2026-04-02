% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function data = All_Torques(params,n)
range_yaw = linspace(-180,180,n);
range_pitch = linspace(-90,90,n);
range_roll = linspace(-180,180,n);
i = 1;

data = cell(1,n^3);

for y = 1:n
    yaw = range_yaw(y);
    for p = 1:n
        pitch = range_pitch(p);
        for r = 1:n
            roll = range_roll(r);
            EA = [yaw,pitch,roll];
            SN = Euler3212C(EA)';
            [ ~, ~, B1_L1, B2_L2, ~, overlapFlag] = ...
                multisphereFT( params.debris.N_spheres, params.servicer.N_spheres, params.N_rvec_km*1000,...
                params.V, eye(3), SN, params.debris.D_COM, params.servicer.S_COM);
            data{i}.EA = EA;
            data{i}.C2 = SN;
            data{i}.B2_L2 = B2_L2;
            data{i}.overlapFlag = overlapFlag;
            i = i+1;
        end
    end
end

end
