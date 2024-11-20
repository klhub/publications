function [q]=q_mult(qa,qb)

% function [q]=q_mult(qa,qb)
% 
% to perform quaternion multiplication
%
% INPUT: -
% qa        first quaternion [4Xm]
% qb        second quaternion [4Xm]
%
% OUTPUT: -
% q        new quaternion [4Xm]    

qa1 = qa(1,:);
qa2 = qa(2,:);
qa3 = qa(3,:);
qa4 = qa(4,:);
qb1 = qb(1,:);
qb2 = qb(2,:);
qb3 = qb(3,:);
qb4 = qb(4,:);

q1 = +qb4.*qa1-qb3.*qa2+qb2.*qa3+qb1.*qa4 ;
q2 = +qb3.*qa1+qb4.*qa2-qb1.*qa3+qb2.*qa4 ;
q3 = -qb2.*qa1+qb1.*qa2+qb4.*qa3+qb3.*qa4 ;
q4 = -qb1.*qa1-qb2.*qa2-qb3.*qa3+qb4.*qa4 ;

q  = [q1;q2;q3;q4];