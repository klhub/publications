function [om_out]=om(w,dt)

%[om_out]=om(w)
%
% for use in discrete quaternion propagation
%
% The inputs are:
%      w = angular velocity [3X1]
%     dt = time step
%
% The outputs are:
% om_out = estimated quaternions [mx4] 

wn  = norm(w);
psi = sin(.5*wn*dt)*w/wn;
Z   = cos(.5*wn*dt)*eye(3)-cross(psi);

om_out = [  Z    psi
           -psi' cos(.5*wn*dt) ];
