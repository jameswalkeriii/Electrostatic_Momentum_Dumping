% Computing all the torque vectors generated between two charged
% spacecraft for all orientations
function data = All_Torques(params,n)
range_yaw = linspace(-180,180,n);
range_pitch = linspace(-90,90,n);
range_roll = linspace(-180,180,n);
i = 1;

for y = 1:n
    yaw = range_yaw(y);
    for p = 1:n
        pitch = range_pitch(p);
        for r = 1:n
            roll = range_roll(r);
            EA = [yaw,pitch,roll];
            C2 = Euler3212C(EA)';
            [ ~, ~, ~, T_L2, ~, overlapFlag] = ...
                multisphereFT( params.debris.spheres, params.servicer.spheres, params.r_km*1000,...
                params.V, eye(3), C2, params.debris.COM, params.servicer.COM);
            data{i}.EA = EA;
            data{i}.C2 = C2;
            B_L2 = C2*T_L2;
            data{i}.B_L = B_L2;
            data{i}.T_L = T_L2;
            data{i}.overlapFlag = overlapFlag;
            i = i+1;
        end
    end
end

end

