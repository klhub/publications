function [mea]=mea_star(q,q_b_s,noise)

% function [mea]=mea_star(q,noise_st)
% 
% to simulate quaternion measurement of star tracker
%
% INPUT: -
% q         current true quaternion, from ECI to body coordinate
% q_b_s     quaternion from body to sensor (true, with misalignment) coordinate
% noise     star tracker noise
%
% OUTPUT: -
% mea       star tracker quaternion measurement    

%q_mea = q_mult(q_b_s,q)+noise;
%q_mea = q_mea/norm(q_mea);

q_mea = q_mult(q_b_s,q);
noise = [noise(1:3); 1]; noise=noise/norm(noise);
q_mea = q_mult(noise,q_mea);
%q_mea = q_mea/norm(q_mea);

mea = q_mea ;
