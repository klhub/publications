function [q_dot]=q_fun(q,w,noise)

q_dot  = .5*omega([w ; 0])*q ;