classdef constants
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
   properties (Constant)
       Kc = 8.9875517923e9;      % coulomb constant [N m^2 C^-2]
       muEarth = 3.986004418e14*1.e-9; % gravitational parameter earth [km^3 s^-2]
       GravConst = 6.67430e-11*1.e-9; % gravitational constant [km^3 kg^-1 s^-2]
       kB = 8.617333262145e-5; % Boltzmann constant [eV/K]
       q0 = 1.602176634e-19; % elementary charge [C]
       me = 9.1093837015e-31; % mass of electron [kg]
       mp = 1.67262192369e-27; % mass of proton [kg]
       eps0 = 8.8541878128e-12; % permittivity of free space [F m^-1]
   end
    
%     methods
%         function obj = untitled(inputArg1,inputArg2)
%             %UNTITLED Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
%         
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
%     end
end

