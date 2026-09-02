function debris_params = load_Debris(debris_shape,ND_i)
% Rotation Matrix for a -90 deg rotation about the Y axis
M_n90deg_Y = [cosd(-90),0,-sind(-90);0,1,0;sind(-90),0,cosd(-90)];

% Rotation Matrix for a 180 deg rotation about the Y axis
M_180deg_Y = [cosd(180),0,-sind(180);0,1,0;sind(180),0,cosd(180)];

% Rotation Matrix for a 180 deg rotation about the Z axis
M_180deg_Z = [cosd(180) , sind(180), 0;-sind(180), cosd(180), 0;0, 0, 1];

% Rotation Matrix for a 90 deg rotation about the Z axis
M_90deg_Z = [cosd(90) , sind(90), 0;-sind(90), cosd(90), 0;0, 0, 1];

% Rotation Matrix for a 90 deg rotation about the X axis
M_90deg_X = [1,0,0;0,cosd(90),-sind(90);0,sind(90),cosd(90)];

% Rotation Matrix for a 90 deg rotation about the X axis
M_180deg_X = [1,0,0;0,cosd(180),-sind(180);0,sind(180),cosd(180)];


switch debris_shape
    case "SSL1300_bus"
        debris_params.shape = "SSL1300_bus";
        debris_params.mass = 2000; % mass [kg];
        debris_params.D_COM = [0, 1.4, 0]'; % [m] SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
        debris_params.D_MI = [2333; 43092; 42725].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan, in body frame
        debris_params.sphs = load('SSL1300_bus.mat');
        
        DI = M_n90deg_Y*M_90deg_X;
               
    case "small_SSL1300_bus_updated"
        % https://handwiki.org/wiki/Engineering%3AIntelsat_706?utm_source=chatgpt.com
        debris_params.shape = "small_SSL1300_bus_updated";
        debris_params.mass = 1450; % mass [kg];
        debris_params.D_COM = [0, -1.1/4, -2.45/4]'; % [m] SSL-1300 ish COM in body frame
        debris_params.D_MI = [1531; 20011; 19473].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan, in body frame
        debris_params.sphs = load('small_SSL1300_bus_updated.mat');
        
        DI = M_n90deg_Y;

    case "big_SSL1300_bus_updated"
        %https://www.researchgate.net/publication/51219764_Navy_Prototype_Optical_Interferometer_observations_of_geosynchronous_satellites
        debris_params.shape = "big_SSL1300_bus_updated";
        debris_params.mass = 2364; % mass [kg];
        debris_params.D_COM = [0, -2.9/6, 0]'; % [m] SSL-1300 ish COM. Distance is from center-front of body (docking location), in body frame
        debris_params.D_MI = [2739; 79706; 79371].*eye(3); % [kg m2] SSL Moment of Inertia From email with Dan, in body frame
        debris_params.sphs = load('big_SSL1300_bus_updated.mat');
        
        DI = M_n90deg_Y;               
        
    case "GOESR_bus_noboom"
        debris_params.shape = "GOESR_bus_noboom";
        debris_params.mass= 2857-25; % [kg] GOES-R series dry mass
        
        debris_params.D_MI = [16104.5-1.26,8216-533.33,12577.1-533.33]'.*eye(3); % in body frame no boom
        debris_params.D_COM = [-87.5,2060,393.7]'*1e-3; % in body frame
        
        debris_params.sphs = load('GOESR_bus_noboom'); % Loading MSM model for GOES-R without boom
        DI = M_n90deg_Y*M_90deg_X;
        
    case "GOESR_bus" %https://www.goes-r.gov/downloads/resources/documents/GOES-RSeriesDataBook.pdf
        debris_params.shape = "GOESR_bus";
        debris_params.mass= 2857; % [kg] GOES-R series dry mass
        
        debris_params.D_MI = [16104.5,8216,12577.1]'.*eye(3); % in body frame
        debris_params.D_COM = [-87.5,2060,393.7]'*1e-3; % in body frame
        
        debris_params.sphs = load('GOESR_bus'); % Loading MSM model for GOES-R without boom
        DI = M_n90deg_Y*M_90deg_X;
        
    case "GOESR_bus_updated"
        debris_params.shape = "GOESR_bus_updated";
        debris_params.mass= 2857; % [kg] GOES-R series dry mass
        
        debris_params.D_MI = [16104.5,8216,12577.1]'.*eye(3); % in body frame
        debris_params.D_COM = [-87.5,393.7,2060]'*1e-3; % in body frame
        
        debris_params.sphs = load('GOESR_bus_updated'); % Loading MSM model for GOES-R without boom
        DI = M_n90deg_Y*M_180deg_X;
        
    case "GOESR_bus_updated_noboom"
        debris_params.shape = "GOESR_bus_updated_noboom";
        debris_params.mass= 2857-25; % [kg] GOES-R series dry mass
        
        debris_params.D_MI = [16104.5-1.26,8216-533.33,12577.1-533.33]'.*eye(3); % in body frame no boom
        debris_params.D_COM = [-87.5,393.7,2060]'*1e-3; % in body frame
        
        debris_params.sphs = load('GOESR_bus_updated_noboom'); % Loading MSM model for GOES-R without boom
        DI = M_n90deg_Y*M_180deg_X;
        
    case "Cylinder"
        debris_params.shape = "Cylinder";
        debris_params.mass = 2000;
        
        debris_params.D_COM = [0,0,0]';
        Ixx = 1/3*debris_params.mass*31.4^2;
        Iyy = Ixx;
        Izz = .5*debris_params.mass*3^2;
        debris_params.D_MI = [Ixx,Iyy,Izz].*eye(3);
        
        cyl = load('31p4m_cylinder','cylinder_31p4m');
        debris_params.sphs.SPHSb = cyl.cylinder_31p4m;
        
        for i = 1:length(debris_params.sphs.SPHSb)
            debris_params.sphs.SPHSb(3,i) = debris_params.sphs.SPHSb(3,i)-31.4/2;
        end
        
        DI = M_90deg_X;
        
    case "Sphere"
        debris_params.shape = "Sphere";
        debris_params.mass = 2000;
        
        debris_params.D_COM = [0,0,0]';
        I = 2/5*debris_params.mass*3^2;
        debris_params.D_MI = [I,I,I].*eye(3);
        
        sph = load('3m_sphere','SPHS_sphere10m');
        debris_params.sphs.SPHSb = sph.SPHS_sphere10m;
        DI = eye(3);
end

for i = 1:length(debris_params.sphs.SPHSb)
    debris_params.D_spheres(1:3,i) = ND_i*DI*debris_params.sphs.SPHSb(1:3,i);
    debris_params.D_spheres(4,i) = debris_params.sphs.SPHSb(4,i);
end

debris_params.D_COM = ND_i*DI*debris_params.D_COM;
debris_params.D_MI = ND_i*DI*debris_params.D_MI*(ND_i*DI)'; % Mapping Momentum of inertia to the new frame


end