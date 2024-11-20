function [x_dot]=states_fun(x,w,noise_q,noise_bd)

q=x(1:4);
b=x(5:7);
q_dot=0.5*xi(q)*w;  % no process noise yet
b_dot=noise_bd;     % bias drift

x_dot=[q_dot;b_dot];
