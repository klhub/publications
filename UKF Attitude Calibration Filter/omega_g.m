function [omega_g]=omega_g(w)

% function [x_dot]=states_true(x,w,n2_noise)
% 
% to propagate true  :  quaternions         , x(1:4)
%                       gyros biases        , x(5:7)
%                       Earth sensor biases , x(8:10)
%

%omega_g   = [  0   w(3) w(2) w(1)  0    0   sign(w(1))     0           0 
%              w(3)  0    0    0   w(2)  0       0      sign(w(2))      0 
%               0    0    0    0    0   w(3)     0          0       sign(w(3)) ] ;

omega_g   = [  0   w(3) w(2) w(1)  0    0  
              w(3)  0    0    0   w(2)  0  
               0    0    0    0    0   w(3) ] ;
