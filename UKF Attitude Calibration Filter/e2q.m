function q=e2q(phi,theta,psi,flag)
%function q=e2q(phi,theta,psi,flag)
%
% This m-file transforms Euler angle to 
% quaternions with the following rotations.
% Theta, phi, and psi are in radians.
%
%  The inputs are:
%    flag = 1 for 1-2-1
%         = 2 for 2-3-2
%         = 3 for 3-1-3
%         = 4 for 1-3-1
%         = 5 for 2-1-2
%         = 6 for 3-2-3
%         = 7 for 1-2-3
%         = 8 for 2-3-1
%         = 9 for 3-1-2
%         = 10 for 1-3-2
%         = 11 for 2-1-3
%         = 12 for 3-2-1

% John L. Crassidis 4/24/95

t1=phi;t2=theta;t3=psi;

if flag==1 | flag==2 | flag==3 | flag==4 | flag==5 | flag ==6
 q4=cos(t2/2).*cos((t1+t3)/2);
end

if flag==1,
 q1=cos(t2/2).*sin((t1+t3)/2);
 q2=sin(t2/2).*cos((t1-t3)/2);
 q3=sin(t2/2).*sin((t1-t3)/2);
elseif flag==2,
 q1=sin(t2/2).*sin((t1-t3)/2);
 q2=cos(t2/2).*sin((t1+t3)/2);
 q3=sin(t2/2).*cos((t1-t3)/2);
elseif flag==3,
 q1=sin(t2/2).*cos((t1-t3)/2);
 q2=sin(t2/2).*sin((t1-t3)/2);
 q3=cos(t2/2).*sin((t1+t3)/2);
elseif flag==4,
 q1=cos(t2/2).*sin((t1+t3)/2);
 q2=sin(t2/2).*sin((t3-t1)/2);
 q3=sin(t2/2).*cos((t3-t1)/2);
elseif flag==5,
 q1=sin(t2/2).*cos((t3-t1)/2);
 q2=cos(t2/2).*sin((t3+t1)/2);
 q3=sin(t2/2).*sin((t3-t1)/2);
elseif flag==6,
 q1=sin(t2/2).*sin((t3-t1)/2);
 q2=sin(t2/2).*cos((t3-t1)/2);
 q3=cos(t2/2).*sin((t3+t1)/2);
end

if flag==7 | flag==8 | flag==9
 q4=cos(t1/2).*cos(t2/2).*cos(t3/2)-sin(t1/2).*sin(t2/2).*sin(t3/2);
end
if flag==10 | flag==11 | flag==12
 q4=cos(t1/2).*cos(t2/2).*cos(t3/2)+sin(t1/2).*sin(t2/2).*sin(t3/2);
end

if flag==7,
 q1=sin(t1/2).*cos(t2/2).*cos(t3/2)+cos(t1/2).*sin(t2/2).*sin(t3/2);
 q2=cos(t1/2).*sin(t2/2).*cos(t3/2)-sin(t1/2).*cos(t2/2).*sin(t3/2);
 q3=cos(t1/2).*cos(t2/2).*sin(t3/2)+sin(t1/2).*sin(t2/2).*cos(t3/2);
elseif flag==8,
 q1=cos(t1/2).*cos(t2/2).*sin(t3/2)+sin(t1/2).*sin(t2/2).*cos(t3/2);
 q2=sin(t1/2).*cos(t2/2).*cos(t3/2)+cos(t1/2).*sin(t2/2).*sin(t3/2);
 q3=cos(t1/2).*sin(t2/2).*cos(t3/2)-sin(t1/2).*cos(t2/2).*sin(t3/2);
elseif flag==9,
 q1=cos(t1/2).*sin(t2/2).*cos(t3/2)-sin(t1/2).*cos(t2/2).*sin(t3/2);
 q2=cos(t1/2).*cos(t2/2).*sin(t3/2)+sin(t1/2).*sin(t2/2).*cos(t3/2);
 q3=sin(t1/2).*cos(t2/2).*cos(t3/2)+cos(t1/2).*sin(t2/2).*sin(t3/2);
elseif flag==10,
 q1=sin(t1/2).*cos(t2/2).*cos(t3/2)-cos(t1/2).*sin(t2/2).*sin(t3/2);
 q2=cos(t1/2).*cos(t2/2).*sin(t3/2)-sin(t1/2).*sin(t2/2).*cos(t3/2);
 q3=cos(t1/2).*sin(t2/2).*cos(t3/2)+sin(t1/2).*cos(t2/2).*sin(t3/2);
elseif flag==11,
 q1=cos(t1/2).*sin(t2/2).*cos(t3/2)+sin(t1/2).*cos(t2/2).*sin(t3/2);
 q2=sin(t1/2).*cos(t2/2).*cos(t3/2)-cos(t1/2).*sin(t2/2).*sin(t3/2);
 q3=cos(t1/2).*cos(t2/2).*sin(t3/2)-sin(t1/2).*sin(t2/2).*cos(t3/2);
elseif flag==12,
 q1=cos(t1/2).*cos(t2/2).*sin(t3/2)-sin(t1/2).*sin(t2/2).*cos(t3/2);
 q2=cos(t1/2).*sin(t2/2).*cos(t3/2)+sin(t1/2).*cos(t2/2).*sin(t3/2);
 q3=sin(t1/2).*cos(t2/2).*cos(t3/2)-cos(t1/2).*sin(t2/2).*sin(t3/2);
end

q=[q1 q2 q3 q4];

