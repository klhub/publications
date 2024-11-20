function [mea]=mea_gyro(w,M,T_g0_b,b,noise)

% function [mea]=mea_gyro(w,C,dd)
% 
% to simulate measurement of skewed gyro quadruplet
%
% INPUT: -
% w     current rotational rate
% C     gyro mounting matrix
% dd    8X1 misalignment, 4X1 scale factor, 4X1 gyro bias
%
% OUTPUT: -
% mea   gyro measurement

mea = inv(eye(3)+M)*T_g0_b*w-b-noise ;
%mea = (eye(3)-M)*T_b_g0*w-b-noise ;
