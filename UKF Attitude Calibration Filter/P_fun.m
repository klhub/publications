function [P_dot]=P_fun(P,F,G,Q)

% This program propagate truth quaternion and generate measurement
% 
% INPUTS: -
% dt    truth propagation time
% 
% OUTPUTS: -
% q     truth quaternion

P_dot = F*P+P*F'+G*Q*G';