function [mea]=mea_star2(q,T_s_b,pr,noise)

% function [mea]=mea_star2(q,T_b_c,pr,noise)
% 
% to simulate unit vector measurement of star tracker
%
% INPUT: -
% q         current true quaternion, from ECI to body coordinate
% T_b_c     tranformation matrix from body to sensor (true, with misalignment) coordinate
% pr        reference vector
% noise     star tracker noise
%
% OUTPUT: -
% mea       star tracker unit vector measurement    

p_mea = T_s_b*q2att_mat(q)*pr+noise;
%p_mea = T_b_s*q2att_mat(q)*pr;
%p_mea = p_mea/norm(p_mea);

mea = p_mea ;
