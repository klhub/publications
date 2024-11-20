function [m]=m_fun(g,w,b)

%  !!! This file is not in used, delete if don't want to use this in future

% function [m]=m_fun(DEL,LAM,UUU,w)
% 
% to propagate true  :  quaternions         , x(1:4)
%                       gyros biases        , x(5:7)
%                       Earth sensor biases , x(8:10)
%

DEL = [ 0 g(3) g(2)
        0  0   g(1)
        0  0    0   ];
LAM = diag(g(4:6));
UUU = diag(g(7:9))*diag(sign(w));

m = DEL+LAM+UUU;
