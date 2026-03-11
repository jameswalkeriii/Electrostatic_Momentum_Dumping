function Xdot = N_RW_EOM(X, I_RW, Iws, Gs0, N,us,L)
sig_BN = X(1:3);
w_BN = X(4:6);

tt = 1;
Om_mat = zeros(N,1);
while tt < N+1
    Om_mat(tt) = X(6+tt);
    tt = tt + 1;
end

Gs = Gs0;
hs = [];
for k = 1:N
    
    gsi = Gs(1:3,k);  
    ws = dot(w_BN,gsi);  
    hsi = Iws*(ws + Om_mat(k));
    hs = [hs;hsi];

end

Bsig = ((1-norm(sig_BN)^2)*eye(3)+2*tild(sig_BN)+2*(sig_BN*sig_BN'));

sig_BN_dot = 1/4*Bsig*w_BN;

w_BN_dot = inv(I_RW)*(-tild(w_BN)*I_RW*w_BN - tild(w_BN)*(Gs*hs)-Gs*us+L);

Om_dot = us/Iws;
% Wheel speed saturation
% for i = 1:length(Om_mat)
%     if Om_mat(i) > 50
%         Om_dot(i) = 0;
%     end
% end

Xdot = [sig_BN_dot;w_BN_dot;Om_dot];
end