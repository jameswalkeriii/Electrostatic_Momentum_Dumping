classdef storage
    properties
        Ttot
        Xtot_servicer
        Xtot_debris
        Htot_servicer
        Htot_debris
        aterrtot_servicer
        werrtot_servicer
        ustot_servicer
        Lrtot_servicer
        L_e_servicer_tot
        L_e_debris_tot
        normBN_servicer
        reftot_servicer
        N_rvec_tot
        flags
        debris_ang_vel
        
    end

   methods
       function obj = storage(Ttot)
           tn = length(Ttot);
           obj.Ttot = Ttot;
           obj.Xtot_servicer = zeros(9,tn);
           obj.Xtot_debris = zeros(6,tn);
           obj.Htot_servicer = zeros(3,tn);
           obj.Htot_debris = zeros(3,tn);
           obj.aterrtot_servicer = zeros(3,tn);
           obj.werrtot_servicer = zeros(3,tn);
           obj.ustot_servicer = zeros(3,tn);
           obj.Lrtot_servicer = zeros(3,tn);
           obj.L_e_servicer_tot = zeros(3,tn);
           obj.L_e_debris_tot = zeros(3,tn);
           obj.normBN_servicer = zeros(1,tn);
           obj.reftot_servicer = zeros(3,tn);
           obj.N_rvec_tot = zeros(3,tn);
           obj.flags = zeros(2,tn);
           obj.debris_ang_vel = zeros(1,tn);
           
       end
      
      function obj = update_storage(obj,tt,params)
         obj.Htot_servicer(:,tt) = params.sim.H;
         obj.Htot_debris(:,tt) = params.sim.H_deb;
         obj.aterrtot_servicer(:,tt) = params.sim.aterr;
         obj.Xtot_servicer(:,tt) = params.sim.X_serv;
         obj.Xtot_debris(:,tt) = params.sim.X_deb;
         obj.werrtot_servicer(:,tt) = params.sim.werr;
         obj.ustot_servicer(:,tt) = params.sim.us;
         obj.Lrtot_servicer(:,tt) = params.sim.Lr;
         obj.L_e_servicer_tot(:,tt) = params.sim.L_e_serv;
         obj.L_e_debris_tot(:,tt) = params.sim.L_e_deb;
         obj.normBN_servicer(:,tt) = norm(params.sim.X_serv(1:3));
         obj.reftot_servicer(:,tt) = params.sim.sig_RN;
         obj.N_rvec_tot(:,tt) = params.sim.N_rvec_m;
         obj.flags(:,tt) = [params.desat_flag; params.sim.mode_code];
         obj.debris_ang_vel(:,tt) = norm(params.sim.X_deb(4:6));
      end
   end
end
