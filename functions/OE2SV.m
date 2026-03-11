function [ r,v ] = OE2SV( a,inclination,argPer,RAAN,theta,enorm)
%orbital elements at initial time
% a = r1; %semi-major axis [kilometers] 
% enorm  = 0; %eccentricity. 
% inclination = 0; %inclination [degrees]
% RAAN = 0; %  right ascension of the ascending node [degrees]
% argPer = 0; % argument of perigee [degrees]
% theta = 0; %True anomaly [degrees]

if nargin == 1
    OE = a;
    a = OE(1);
    inclination = OE(2);
    argPer = OE(3);
    RAAN = OE(4);
    theta = OE(5);
    enorm = OE(6);
end

mu = constants.muEarth;

%converting orbital parameters to ECI cartsian coordinates 
p = a*(1-enorm^2); %intermediate variable
q = p/(1+enorm*cosd(theta));%intermediate variable

% Creating r vector in pqw coordinates
R_pqw(1,1) = q*cosd(theta);
R_pqw(2,1) = q*sind(theta);
R_pqw(3,1) = 0;
    
% Creating v vector in pqw coordinates    
V_pqw(1,1) = -(mu/p)^.5*sind(theta);
V_pqw(2,1) = ((mu/p)^.5)*(enorm + cosd(theta));
V_pqw(3,1) =   0;

%Solving for 313 rotation matrices
R1_i = [1 0 0; 0 cosd(inclination) sind(inclination); 0 -sind(inclination) cosd(inclination)];
R3_Om = [cosd(RAAN) sind(RAAN) 0; -sind(RAAN) cosd(RAAN) 0; 0 0 1];
R3_om = [cosd(argPer) sind(argPer) 0; -sind(argPer) cosd(argPer) 0; 0 0 1];

support_var = R3_om*R1_i*R3_Om; %Intermediate variable

r = support_var'*R_pqw; %Radius r [km] in ECI Cartesian
v = support_var'*V_pqw; %Velocity v [km/s] in ECI Cartesian


end

