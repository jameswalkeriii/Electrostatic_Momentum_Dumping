classdef storage
    properties
        Ttot
        Xtot
        Htot
        aterrtot
        werrtot
        ustot
        Lrtot
        L_e_tot
        normBN
        reftot
        rtot
        flags
        
    end

   methods
       function obj = storage(Ttot)
           tn = length(Ttot);
           obj.Ttot = Ttot;
           obj.Htot = zeros(3,tn);
           obj.Xtot = zeros(9,tn);
           obj.aterrtot = zeros(3,tn);
           obj.werrtot = zeros(3,tn);
           obj.ustot = zeros(3,tn);
           obj.Lrtot = zeros(3,tn);
           obj.L_e_tot = zeros(3,tn);
           obj.normBN = zeros(1,tn);
           obj.reftot = zeros(3,tn);
           obj.rtot=zeros(3,tn);
           obj.flags = zeros(2,tn);
           
       end
      
      function obj = update_storage(obj,tt,params)
         obj.Htot(:,tt) = params.sim.H;
         obj.aterrtot(:,tt) = params.sim.aterr;
         obj.Xtot(:,tt) = params.sim.X_servicer;
         obj.werrtot(:,tt) = params.sim.werr;
         obj.ustot(:,tt) = params.sim.us;
         obj.Lrtot(:,tt) = params.sim.Lr;
         obj.L_e_tot(:,tt) = params.sim.L_e;
         obj.normBN(:,tt) = norm(params.sim.X_servicer(1:3));
         obj.reftot(:,tt) = params.sim.sig_RN;
         obj.rtot(:,tt) = params.sim.r_m;
         obj.flags(:,tt) = [params.desat_flag; params.sim.mode_code];
      end
   end
end
