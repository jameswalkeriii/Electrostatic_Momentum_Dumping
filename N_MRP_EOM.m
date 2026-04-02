function Xdot = N_MRP_EOM(X,I,L)
sig_BN = X(1:3);
w_BN = X(4:6);

sig_BN_dot = (1/4)*((1-norm(sig_BN)^2)*eye(3) + 2*tild(sig_BN) + 2*(sig_BN*sig_BN'))*w_BN;
w_BN_dot = inv(I)*(-tild(w_BN)*I*w_BN+L);

Xdot = [sig_BN_dot;w_BN_dot];
end