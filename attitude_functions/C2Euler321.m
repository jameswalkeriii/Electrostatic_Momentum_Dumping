function q = C2Euler321(C)

% C2Euler321
%
%   Q = C2Euler321(C) returns the 3x1 vector of 3-2-1 Euler angles
%   [yaw; pitch; roll] in radians corresponding to the direction cosine
%   matrix C used by Euler3212C.
%
%   At the 3-2-1 singularity (pitch = +/- pi/2), yaw and roll are not
%   independently observable. This function sets yaw = 0 and returns the
%   equivalent roll angle.
%

pitch = asin(-C(1,3));
cp = cos(pitch);
tol = 1e-12;

if abs(cp) > tol
    yaw = atan2(C(1,2), C(1,1));
    roll = atan2(C(2,3), C(3,3));
else
    yaw = 0;

    if C(1,3) <= -1 + tol
        % pitch = +pi/2, only (roll - yaw) is observable
        roll = atan2(C(2,1), C(2,2));
    else
        % pitch = -pi/2, only (roll + yaw) is observable
        roll = atan2(-C(2,1), C(2,2));
    end
end

q = [yaw; pitch; roll];
