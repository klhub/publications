function A = q2att_mat(qua)

q = [ qua(1) ; qua(2) ; qua(3) ; qua(4) ] ;

A = (q(4)^2-q(1:3)'*q(1:3))*eye(3)+2*q(1:3)*q(1:3)'-2*q(4)*cross(q(1:3)) ;
