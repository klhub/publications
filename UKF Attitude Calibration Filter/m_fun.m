function [m]=m_fun(g,w)

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
%UUU = diag(g(7:9))*diag(sign(w));
%UUU = [g(7)*sign(w(1))     0               0 
%           0           g(8)*sign(w(2))     0
%           0               0           g(9)*sign(w(3)) ];

%m = DEL+LAM+UUU;
m = DEL+LAM;
