function [q_diff,w_diff,b_diff,g_diff,s1_diff,p_diff]=err(dts,tf,q,qe,w,we,b,be,ge,s1e,pe);
% [q_diff,w_diff,b_diff,g_diff,s1_diff]=diff(dts,tf,q,qe,w,we,b,be,ge,s1e);

% function [mea]=mea_star(q,noise_st)
% 
% Kalman Filtering for Spacecraft System Alignment Calibration, Mark E. Pittelkau, JGCD Nov.-Dec. '01
%
% INPUT: -
% dts       sampling interval
%
% OUTPUT: -
% qe        estimated quaternion

% constants
parameters;             % loading the sensors misalignments, noise properties ...

% pre-allocate space
q_diff  = zeros(4,tf/dts);      
w_diff  = zeros(3,tf/dts);      
b_diff  = zeros(3,tf/dts); 

g_diff  = zeros(6,tf/dts);
s1_diff = zeros(3,tf/dts);
p_diff  = zeros(3,tf/dts);

% difference
g  = [ kron(ones(1,tf/dts),g(1))
       kron(ones(1,tf/dts),g(2))
       kron(ones(1,tf/dts),g(3))
       kron(ones(1,tf/dts),g(4))
       kron(ones(1,tf/dts),g(5))
       kron(ones(1,tf/dts),g(6)) ];
%       kron(ones(1,tf/dts),g(7))
%       kron(ones(1,tf/dts),g(8))
%       kron(ones(1,tf/dts),g(9)) ] ;
s1 = [ kron(ones(1,tf/dts),s1(1))  
       kron(ones(1,tf/dts),s1(2))
       kron(ones(1,tf/dts),s1(3)) ];        % QUATERNION ST
s1 = [ kron(ones(1,tf/dts),del_st1(1))  
       kron(ones(1,tf/dts),del_st1(2))
       kron(ones(1,tf/dts),del_st1(3)) ];   % VECTOR ST
pp = [ kron(ones(1,tf/dts),del_p(1))  
       kron(ones(1,tf/dts),del_p(2))
       kron(ones(1,tf/dts),del_p(3)) ];   % VECTOR payload
g_diff  = ge-g;
s1_diff = s1e-s1;
p_diff  = pe-pp;

wait = waitbar(0,'Difference running ...');
for i=1:tf/dts-1
    if i==1
        w_diff(:,i) = we(:,i)-w(:,2*i-1);
        b_diff(:,i) = be(:,i)-b(:,2*i-1);
    end
    q_diff(:,i) = delta_q(qe(:,i),q(:,2*i-1));
    if i>1
        w_diff(:,i) = we(:,i)-w(:,2*i-3);
        b_diff(:,i) = be(:,i)+b(:,2*i-3);
        %q_diff(:,i) = delta_q(qe(:,i),q(:,2*i-1));
    end

    %a=delta_q(qe(:,i),q(:,2*i-1)); [aa,bb,cc]=q2e(a',7); q_diff(1:3,i)=[aa;bb;cc];     % "exact" error, but almost identical results
    %if i>1                                                                             % if use this, remember to take out the 2*
    %    a=delta_q(qe(:,i),q(:,2*i-3)); [aa,bb,cc]=q2e(a',7); q_diff(1:3,i)=[aa;bb;cc]; % in plotting the attitude error in deg
    %end
    
    waitbar(i/(tf/dts),wait)
end
close(wait)

% Just to make the plots slightly nicer, avoid those 
q_diff(:,1)=q_diff(:,5);   q_diff(:,2)=q_diff(:,5);   q_diff(:,3)=q_diff(:,5);   q_diff(:,4)=q_diff(:,5);   q_diff(:,tf/dts)=q_diff(:,tf/dts-1);
w_diff(:,1)=w_diff(:,5);   w_diff(:,2)=w_diff(:,5);   w_diff(:,3)=w_diff(:,5);   w_diff(:,4)=w_diff(:,5);   w_diff(:,tf/dts)=w_diff(:,tf/dts-1); 
b_diff(:,1)=b_diff(:,5);   b_diff(:,2)=b_diff(:,5);   b_diff(:,3)=b_diff(:,5);   b_diff(:,4)=b_diff(:,5);   b_diff(:,tf/dts)=b_diff(:,tf/dts-1); 
g_diff(:,1)=g_diff(:,3);   g_diff(:,2)=g_diff(:,5);   g_diff(:,3)=g_diff(:,5);   g_diff(:,4)=g_diff(:,5);   g_diff(:,tf/dts)=g_diff(:,tf/dts-1); 
s1_diff(:,1)=s1_diff(:,5); s1_diff(:,2)=s1_diff(:,5); s1_diff(:,1)=s1_diff(:,5); s1_diff(:,3)=s1_diff(:,5); s1_diff(:,4)=s1_diff(:,5); s1_diff(:,1)=s1_diff(:,3); s1_diff(:,tf/dts)=s1_diff(:,tf/dts-1); 
