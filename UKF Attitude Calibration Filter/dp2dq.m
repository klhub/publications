function [dq]=dp2dq(dp,fac,par)

%[om_out]=om(w)
%
% for use in discrete quaternion propagation
%
% The inputs are:
%     dp = error generalized Rodrigues parameters
%    fac = f, scale factor
%    par = a, a parameter from 0 to 1
%
% The outputs are:
%     dq = error quaternion 

dpsq  = dp'*dp;   % dp square
qerr4 = (fac*sqrt(fac^2+(1-par^2)*dpsq)-par*dpsq)/(fac^2+dpsq);
dq    = [ (par+qerr4)/fac*dp ; qerr4 ];
