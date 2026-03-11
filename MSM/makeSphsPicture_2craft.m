function makeSphsPicture_2craft( SPHS1, SPHS2, dR1, dR2, V, C1, C2, C)

% Function to plot two MSM models with proper sphere radii and coloring to
% represent the charge on each sphere. 
% Inputs:
% SPHS1 = [x1 x2 ...
%          y1 y2 ...
%          z1 z2 ...
%          R1 R2 ...]
%        (positions & radii of spheres composing body 1)
% SPHS2 = (same for body 2)
% C1/2 = DCM of body 1/2 rotation
% dR1/dR2 = [xr;yr;zr] (translational posiition of each body)
% V = [V1 V2] (voltages of bodies 1 and 2)

if nargin <6
    C1 = eye(3);
    C2 = eye(3);
end
    SPHS1(1:3,:) = C1*SPHS1(1:3,:); %Rotate spheres by DCMs
    x = SPHS1(1,:); 
    y = SPHS1(2,:); 
    z = SPHS1(3,:); 
        
    R = SPHS1(4,:);
    
    Ns = numel(x);

    % build Sc matrix
    Sc = zeros(Ns);
    for i=1:Ns
        Sc(i,i)=1/SPHS1(4,i);
        for j=i+1:Ns
            Sc(i,j)=1/norm(SPHS1(1:3,j)-SPHS1(1:3,i));
            Sc(j,i)=Sc(i,j);
        end
    end
    
    % build Vs and solve for qs
    k = constants.Kc;
    Vs = V(1) *ones(Ns,1);
    qc = (k*Sc)\(Vs);
        
    [xs,ys,zs] = sphere(15);
%     h = figure;
    hold on

    for i = 1:Ns
        xp = R(i).*xs + x(i)*ones(size(xs))+dR1(1); 
        yp = R(i).*ys + y(i)*ones(size(xs))+dR1(2);
        zp = R(i).*zs + z(i)*ones(size(xs))+dR1(3);
        if nargin < 8
            color = qc(i)*ones(size(xp))*10^9;
        else
            color = C*ones(size(xp));
        end
        surf(xp, yp, zp, color,'FaceAlpha',0.75, 'EdgeColor','none')
    end
    
%% Second craft:
    % rotate positions of spheres according to DCMs
    SPHS2(1:3,:) = C2*SPHS2(1:3,:);
    x = SPHS2(1,:); 
    y = SPHS2(2,:); 
    z = SPHS2(3,:); 
    R = SPHS2(4,:);
    
    Ns = numel(x);

    % build Sc matrix
    Sc = zeros(Ns);
    for i=1:Ns
        Sc(i,i)=1/SPHS2(4,i);
        for j=i+1:Ns
            Sc(i,j)=1/norm(SPHS2(1:3,j)-SPHS2(1:3,i));
            Sc(j,i)=Sc(i,j);
        end
    end
    
    % build Vs and solve for qs
    k = constants.Kc;
    Vs = V(2) *ones(Ns,1);
    qc = (k*Sc)\(Vs);
        
    [xs,ys,zs] = sphere(15);
    hold on
    
    for i = 1:Ns
        xp = R(i).*xs + x(i)*ones(size(xs))+dR2(1); 
        yp = R(i).*ys + y(i)*ones(size(xs))+dR2(2);
        zp = R(i).*zs + z(i)*ones(size(xs))+dR2(3);

        if nargin < 8
            color = qc(i)*ones(size(xp))*10^9;
        else
            color = C*ones(size(xp));
        end
        surf(xp, yp, zp, color,'FaceAlpha',0.75, 'EdgeColor','none')
    end

    axis equal
    grid on
    grid on
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')

end

