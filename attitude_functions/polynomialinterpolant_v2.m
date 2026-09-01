function refTraj = polynomialinterpolant(Xi,Yi,Zi,VXi,VYi,VZi,thetai, dthetai, alpha,beta,gamma,refTraj)
% Computes the time vectors of Y and Z given the initial position and
% velocity, the last n-2 polynomial coefficients of Y and Z, and the
% reference X trajectory (refTraj). The polynomials are defined as
%
% Y(s) = alpha2*X^2 + alpha3*X^3 + alpha4*X^4 + ...
%
% so that the initial condition (Y0, VY0) and the final state (0,0) are 
% imposed by 2 linear equations in terms of the coefficients (see notes) 

% Extract values from refTraj
X     = refTraj.refTraj(1,:);
X0    = max(X);
xi    = Xi/X0;
dX    = refTraj.refTraj(4,:);

% Polynomial terms, determined by imposing Y(0) and Y(1) and their
% derivatives
alpha = [0, 0, alpha];
beta  = [0, 0, beta];
gamma = [0, 0, gamma];
ny    = length(alpha);
nz    = length(beta);
ntheta= length(gamma);

if 1 %Xi >=X0-0.1 % Disconnect derivative tracking when we are close to the origin
    if VXi == 0
        alpha(1:2) = ([xi^2, xi^3; 2*xi, 3*xi^2]\[Yi - sum(alpha(3:end).*xi.^(4:(length(alpha)+1)));...
            0 * X0 - sum((4:(length(alpha)+1)) .* alpha(3:end) .* xi.^(3:(length(alpha))))])';
        beta(1:2)  = ([xi^2, xi^3; 2*xi, 3*xi^2]\[Zi - sum(beta(3:end) .*xi.^(4:(length(beta)+1)));...
            0 * X0 - sum((4:(length(beta)+1)) .* beta(3:end) .* xi.^(3:(length(beta))))])';
        gamma(1:2) = ([xi^2, xi^3; 2*xi, 3*xi^2]\[thetai - sum(gamma(3:end) .*xi.^(4:(length(gamma)+1)));...
            0 * X0 - sum((4:(length(gamma)+1)) .* gamma(3:end) .* xi.^(3:(length(gamma))))])';
    else
        alpha(1:2) = ([xi^2, xi^3; 2*xi, 3*xi^2]\[Yi - sum(alpha(3:end).*xi.^(4:(length(alpha)+1)));...
            VYi/VXi * X0 - sum((4:(length(alpha)+1)) .* alpha(3:end) .* xi.^(3:(length(alpha))))])';
        beta(1:2)  = ([xi^2, xi^3; 2*xi, 3*xi^2]\[Zi - sum(beta(3:end) .*xi.^(4:(length(beta)+1)));...
            VZi/VXi * X0 - sum((4:(length(beta)+1)) .* beta(3:end) .* xi.^(3:(length(beta))))])';
        gamma(1:2) = ([xi^2, xi^3; 2*xi, 3*xi^2]\[thetai - sum(gamma(3:end) .*xi.^(4:(length(gamma)+1)));...
            dthetai/VXi * X0 - sum((4:(length(gamma)+1)) .* gamma(3:end) .* xi.^(3:(length(gamma))))])';
    end
else
    alpha(1:2) = ([xi^2, xi^3; 1, 1]\[Yi - sum(alpha(3:end).*xi.^(4:(length(alpha)+1))); - sum(alpha(3:end))])';
    beta(1:2)  = ([xi^2, xi^3; 1, 1]\[Zi - sum(beta(3:end) .*xi.^(4:(length(beta)+1))); - sum(beta(3:end))])';
    gamma(1:2) = ([xi^2, xi^3; 1, 1]\[thetai - sum(gamma(3:end) .*xi.^(4:(length(gamma)+1))); - sum(gamma(3:end))])';
end

% Polynomial vectors
Y      = zeros(1,length(X));
Z      = zeros(1,length(X));
theta  = zeros(1,length(X));
for i = 1:ny
    Y = Y + alpha(i) * (X/X0).^(i+1);
end
for i = 1:nz
    Z = Z + beta(i) * (X/X0).^(i+1);
end
for i = 1:ntheta
    theta = theta + gamma(i) * (X/X0).^(i+1);
end

% Derivatives
dY      = zeros(1,length(X));
dZ      = zeros(1,length(X));
dtheta  = zeros(1,length(X));
for i = 1:ny
    dY = dY + (i+1) * alpha(i) * (X/X0).^i .* dX/X0;
end
for i = 1:nz
    dZ = dZ + (i+1) * beta(i)  * (X/X0).^i .* dX/X0;
end
for i = 1:ntheta
    dtheta = dtheta + (i+1) * gamma(i)  * (X/X0).^i .* dX/X0;
end

% Store reftraj
refTraj.refTraj(2,:) = Y;
refTraj.refTraj(3,:) = Z;
refTraj.refTraj(5,:) = dY;
refTraj.refTraj(6,:) = dZ;
refTraj.refTraj(7,:) = theta;
refTraj.refTraj(8,:) = dtheta;

% figure
% subplot(2,2,1)
% plot(X,Y)
% hold on
% plot(Xi,Yi,'r+')
% grid on
% subplot(2,2,2)
% plot(X,Z)
% hold on
% plot(Xi,Zi,'r+')
% subplot(2,2,3)
% plot(X,dY)
% hold on
% plot(Xi,VYi,'r+')
% grid on
% subplot(2,2,4)
% plot(X,dZ)
% hold on
% plot(Xi,VZi,'r+')
% grid on
% pause(1)
