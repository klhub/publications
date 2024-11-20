function [mea]=mea_payload(q,T_b_c,pr,noise)

% function [mea]=mea_star(q,noise_st)
% 
% to simulate quaternion measurement of star tracker
%
% INPUT: -
% q         current true quaternion, from ECI to body coordinate
% q_b_s     quaternion from body to sensor (true, with misalignment) coordinate
% noise_st  star tracker noise
%
% OUTPUT: -
% mea       star tracker quaternion measurement    

p_mea = T_b_c*q2att_mat(q)*pr+noise; 
%p_mea = T_b_c*q2att_mat(q)*pr; 
%p_mea = p_mea/norm(p_mea);

mea = p_mea ;
