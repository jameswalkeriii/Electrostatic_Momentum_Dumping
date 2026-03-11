function [ F1, F2, L1, L2, qs, overlapFlag] = multisphereFT( SPHS1, SPHS2, r, V, C1, C2, COM1, COM2)
% Find force and torque between two MSM bodies
% Inputs:
% SPHS1 = [x1 x2 ...
%          y1 y2 ...
%          z1 z2 ...
%          R1 R2 ...]
%        (positions & radii of spheres composing body 1)
% SPHS2 = (same for body 2)
% C1/2 = DCM of body 1/2 rotation
% r = [xr;yr;zr] (position of body 2 origin relative to 1)
% V = [V1 V2] (voltages of bodies 1 and 2)
% COM1/2 = Center of mass to origin offset
%
% Outputs:
% F1, F2, L1, L2 = forces and torques (3x1 vectors) on body 1 and 2
% qs = vector of charges on each sphere

% assume zero rotations if not specified
if nargin < 5
    C1 = eye(3);
    C2 = eye(3);
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
elseif nargin < 6
    C2 = eye(3);
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
elseif nargin < 7
    COM1 = [ 0 0 0 ]';   
    COM2 = [ 0 0 0 ]';
end

% Coulomb's constant [Nm^2/C^2]
k = constants.Kc;

% number of spheres in body 1, body 2, total
n1 = size(SPHS1,2);
n2 = size(SPHS2,2);
n = n1 + n2; 

% rotate positions of spheres according to DCMs
SPHS1t = SPHS1;
SPHS1t(1:3,:) = C1*(SPHS1(1:3,:)) + r*ones(1,n1); 
SPHS2(1:3,:) = C2*(SPHS2(1:3,:)); %

% build matrix with all spheres
SPHS = [SPHS1t SPHS2];

% Find Cinv matrix
% needs to evaluate every sphere to find each charge
sph1 = repmat(SPHS(1:3,:),1,size(SPHS,2));
sph2 = repelem(SPHS(1:3,:),1,size(SPHS,2));

relPos = sum((sph1-sph2).^2,1).^0.5; %vecnorm, but 1/10 the time. Finds relative positions between each pair of spheres in the two bodies

relPos = reshape(relPos, n, n); %distance between each sphere: wil fill off-diagonal elements
Cinv = relPos + SPHS(4,:).*eye(n);

% check if any spheres are overlapping
% lapCheck = (Cinv - repmat(1.1*SPHS(4,:), size(SPHS,2), 1));
% check if distance between sphere1 and 2 is larger than r1 && r2
lapCheck = relPos(n1+1:end, 1:n1) - 1.1*(SPHS1(4,:)+SPHS2(4,:)');

overlapFlag = 0;
% if any(triu(lapCheck,1) < 0, 'all')
if any(lapCheck < 0, 'all')

    overlapFlag = 1;
%     disp('Overlapping spheres')
%     figure; scatter3(SPHS1t(1,:),SPHS1t(2,:),SPHS1t(3,:)); hold on; scatter3(SPHS2(1,:),SPHS2(2,:),SPHS2(3,:),'filled');axis equal
% figure; makeSphsPicture_2craft( SPHS1t, SPHS2, [0 0 0], [0 0 0], V(1), eye(3), eye(3))
end


% figure; makeSphsPicture_2craft( SPHS1t, SPHS2, [0 0 0], [0 0 0], V(1), eye(3), eye(3)); view(3)

Cinv = Cinv.^-1;

% Find charge on each sphere
Vs = [V(1)*ones(n1,1); V(2)*ones(n2,1)];
qs = (k*Cinv)\Vs;

% seperate charges into body 1 & 2
qs1 = qs(1:n1);
qs2 = qs(n1+1:end);

% Find force
% F = k*q1*q2*r/r^3

%compute distances 
sph1 = repmat(SPHS1t(1:3,:), 1, size(SPHS2,2));
sph2 = repelem(SPHS2(1:3,:), 1, size(SPHS1t,2));

relPos = sum((sph1-sph2).^2,1).^0.5;
relVec = sph1-sph2;

Qs1 = repmat(qs1',1,length(qs2));
Qs2 = repelem(qs2',1,length(qs1));

Feach = k*Qs1.*Qs2./relPos.^3;
Feach = Feach.*relVec; %add direction

F1 = sum(Feach,2); 
F2 = -F1; % Force should be equal but opposite

% Find torque
% L = r X f
% compute distances RELATIVE TO CENTER OF MASSES 
sph1 = repmat(SPHS1t(1:3,:) + COM1, 1, size(SPHS2,2));
sph22 = repelem(SPHS2(1:3,:) + COM2, 1, size(SPHS1t,2));

L1 = sum(cross(sph1,Feach),2); %matches L1 to eps
L2 = -sum(cross(sph22,Feach),2); % matches L2 to eps-ish

if norm(L2) > 1 || norm(F2) > 1
%     disp('Oddly high forces/torques')
end

% if norm(abs(L1)-abs(L2)) < eps() && nnz(V)>1
%     disp('Torques are equal')
% end


end



