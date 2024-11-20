function [qe,we,be,ge,s1e,P_cov]=ekf(dts,tf,y)

% function [mea]=mea_star(q,noise_st)
% 
% Kalman Filtering for Spacecraft System Alignment Calibration, Mark E. Pittelkau, JGCD Nov.-Dec. '01
%
% INPUT: -
% dts       sampling interval
% tf        final sampling time
% y         measurements, 3X1 gyros, 4X1 star sensor, 3X1 payload
%
% OUTPUT: -
% qe        estimated quaternion
% we        estimated spacecraft rotational rate
% be        estimated gyro bias
% ge        estimated gyro misalignments, scale factor error
% s1e       estimated star tracker 1 misalignment
% P_cov     output covariance

% constants
T=dts;                  % sampling interval
parameters;             % loading the sensors misalignments, noise properties ...

% pre-allocate space
qe_his  = zeros(4,tf/dts);      % to store all the estimated quaternions, qe
we_his  = zeros(3,tf/dts);      % estimated rotational rate
be_his  = zeros(3,tf/dts);      % estimated gyro bias
ge_his  = zeros(9,tf/dts);      % estimated gyro misalignments, scale factor error
s1e_his = zeros(3,tf/dts);      % estimated star tracker 1 misalignment
P_cov   = zeros(18,tf/dts);     % covariance

% initialization
qe=[0;0;0;1]; be=[0;0;0]*pi/180/3600; %be=bg;    %%%%% !!!!!
ge=[0;0;0;0;0;0;0;0;0]; s1e=[0;0;0]; %ge=g; s1e=s1; %%%%%!!!!!
%we=y(1:3,1)+be;
%Pa=(0.5*pi/180)^2; Pb=(1e-6)^2; Pk=(1500/1e6)^2; Pm=(500/60^2*pi/180)^2; Ps=(60/60^2*pi/180)^2; Pp=(60/60^2*pi/180)^2;
Pa=(0.5*pi/180)^2; Pb=(.5*pi/180/3600)^2; Pk=(500/1e6)^2; Pm=(500/60^2*pi/180)^2; Ps=(50/60^2*pi/180)^2; Pp=(60/60^2*pi/180)^2;
Pa=(5*pi/180)^2; Pb=(.5*pi/180/3600)^2; Pk=(500/1e6)^2; Pm=(500/60^2*pi/180)^2; Ps=(50/60^2*pi/180)^2; Pp=(60/60^2*pi/180)^2;
P = diag([Pa Pa Pa Pb Pb Pb Pm Pm Pm Pk Pk Pk Pk Pk Pk Ps Ps Ps]);
Q = Qc(1:18,1:18);
%R = [SIG_s1m zeros(3); zeros(3) SIG_p];    % QUATERNION output startracker
R = [SIG_s1pm zeros(3); zeros(3) SIG_p];    % VECTOR output startracker

%Pa=0.5*pi/180; Pb=.5*pi/180/3600; Pk=500/1e6; Pm=500/60^2*pi/180; Ps=50/60^2*pi/180; Pp=50/60^2*pi/180;
%P_temp = [Pa;Pa;Pa;Pb;Pb;Pb;Pk;Pm;Pm;Pm;Pk;Pm;Pm;Pm;Pk;Ps;Ps;Ps];
%P=P_temp*P_temp';
P_cov(:,1)=diag(P);

wait = waitbar(0,'Extended Kalman Filter running ...');
for i=1:tf/dts-1
    wgm = y(1:3,i);
    Me=m_fun(ge,wgm);
    we=T_g0_b*[wgm + (eye(3)+Me)*be + omega_g(wgm)*ge];
    % states propagation
    %w=norm(we); co=cos(.5*w*dts); si=sin(.5*w*dts);
    %qw=[we/w;co];
    %qe=omega(qw)*qe;        % propagated quaternions, other states takes the previously updated values
    %qe=expm(.5*dts*omega(we))*qe;       %%%%% !!!!!  % this one works good
    qe=om(we,dts)*qe;       %%%%% !!!!!  % this one works good
    
    %x=[qe;be];
    %F1 = dts*states_fun(x      , we, zeros(4,1), zeros(3,1));
    %F2 = dts*states_fun(x+.5*F1, we, zeros(4,1), zeros(3,1));
    %F3 = dts*states_fun(x+.5*F2, we, zeros(4,1), zeros(3,1));
    %F4 = dts*states_fun(x+F3   , we, zeros(4,1), zeros(3,1));
    %x  = x+(F1+2*F2+2*F3+F4)/6 ;
    %qe = x(1:4);                             % current true quaternions
    
    % covariance propagation
    %Og  = omega_g(wgm+be);
    Og  = omega_g(wgm);    
    G   = [ .5*T_g0_b*(eye(3)+Me) zeros(3,3) zeros(3,9) zeros(3,3)
            zeros(3,3)            eye(3,3)   zeros(3,9) zeros(3,3)
            zeros(9,3)            zeros(9,3) eye(9,9)   zeros(9,3)
            zeros(3,3)            zeros(3,3) zeros(3,9) eye(3,3)   ];
    phi = [ eye(3,3)   .5*T*eye(3)*T_g0_b .5*T*eye(3)*T_g0_b*Og zeros(3,3)  % approximation, but
            zeros(3,3) eye(3,3)           zeros(3,9)            zeros(3,3)  % works great!
            zeros(9,3) zeros(9,3)         eye(9,9)              zeros(9,3)
            zeros(3,3) zeros(3,3)         zeros(3,9)            eye(3,3)   ];

    %tt  = norm(we*T);
    %aa  = cos(tt)*eye(3)+(1-cos(tt))/tt^2*(we*T)*(we*T)'-sin(tt)/tt*cross(we*T);
    %bb  = T*(sin(tt)/tt*eye(3)+(tt-sin(tt))/tt^3*(we*T)*(we*T)'-(1-cos(tt))/tt^2*cross(we*T));
    %phi = [ aa         .5*bb*T_g0_b .5*T*eye(3)*T_g0_b*Og zeros(3,3)    % more "accurate" but performance
    %        zeros(3,3) eye(3,3)     zeros(3,9)            zeros(3,3)    % about the same
    %        zeros(9,3) zeros(9,3)   eye(9,9)              zeros(9,3)
    %        zeros(3,3) zeros(3,3)   zeros(3,9)            eye(3,3)   ];

    %RR  = [ T_g0_b     zeros(3,3) zeros(3,9) zeros(3,3)        % this blocks would decrease the P but not more accurate estimate
    %        zeros(3,3) eye(3,3)   zeros(3,9) zeros(3,3)
    %        zeros(9,3) zeros(9,3) eye(9,9)   zeros(9,3)
    %        zeros(3,3) zeros(3,3) zeros(3,9) eye(3,3)   ];
    %QQ  = [ T/4*SIG_a+T^3/12*SIG_r+T^3/12*Og*SIG_g*Og' T^2/4*SIG_r T^2/4*Og*SIG_g zeros(3,3)
    %        T^2/4*SIG_r                                T*SIG_r     zeros(3,9)     zeros(3,3)
    %        T^2/4*SIG_g*Og'                            zeros(9,3)  T*SIG_g        zeros(9,3)
    %        zeros(3,3)                                 zeros(3,3)  zeros(3,9)     T*SIG_s1   ];
    %Q   = RR*QQ*RR' ;
    P   = phi*P*phi'+G*Q*G';
    
    %F   = [ -cross(we)  .5*T_g0_b  .5*T_g0_b*Og zeros(3,3)  
    %         zeros(3,3) zeros(3,3) zeros(3,9)   zeros(3,3) 
    %         zeros(9,3) zeros(9,3) zeros(9,9)   zeros(9,3) 
    %         zeros(3,3) zeros(3,3) zeros(3,9)   zeros(3,3) ];
    %F1 = dts*P_fun(P      , F, G, Q);
    %F2 = dts*P_fun(P+.5*F1, F, G, Q);
    %F3 = dts*P_fun(P+.5*F2, F, G, Q);
    %F4 = dts*P_fun(P+F3   , F, G, Q);
    %P  = P+(F1+2*F2+2*F3+F4)/6 ;    
    
    % propagated measurement
    %q_s01_s1e=[.5*s1e;1]; q_s01_s1e=q_s01_s1e/norm(q_s01_s1e); % quaternion
    %q_b_s1e=q_mult(q_s01_s1e,q_b_s01); T_b_s1e=q2att_mat(q_b_s1e); % quaternion
    %star1_h    = mea_star(qe,q_b_s1e,zeros(4,1));  % quaternion
    T_s01_s1e    = eye(3)-cross(s1e);  % from nominal sensor mounting to true sensor mounting  %%%% normalize? % VECTOR
    T_b_s1e      = T_s01_s1e*T_b_s01;  % from assumed (nominal) sensor mounting to actual moutning    % VECTOR
    star1_h    = mea_star2(qe,T_b_s1e,pr_st1,zeros(3,1));   % vector
    payload_h  = mea_payload(qe,T_b_c,pr,zeros(3,1));
    %payload_z  = y(8:10,i+1)-payload_h; % payload_z=payload_z/norm(payload_z); % don't normalize    
    payload_z  = y(7:9,i+1)-payload_h; % payload_z=payload_z/norm(payload_z); % don't normalize    % VECTOR
    %star1_z    = xi(y(4:7,i))'*star1_h;        % doesn't work that well
    %star1_z    = delta_q(star1_h,y(4:7,i+1));     % QUATERNION
    star1_z    = y(4:6,i+1)-star1_h;
    z          = [star1_z(1:3);payload_z];  % I think I miscoded it as 2*star1_z(1:3) before
%    z          = [star1_z(1:3)];                        % without payload
    Pr = q2att_mat(qe)*pr; PrX = cross(Pr);             % beware of notation confusion, Pr~=pr
    Pr_st1 = q2att_mat(qe)*pr_st1; Pr_st1X = cross(Pr_st1);             % VECTOR 
    Pc1    = T_b_s1e*q2att_mat(qe)*pr_st1; Pc1X = cross(Pc1);   % VECTOR
    %H  = [ T_b_s1e     zeros(3,12) .5*eye(3)
    %       2*T_b_c*PrX zeros(3,12) zeros(3,3) ];        % sensitivity matrix % QUATERNION
    H  = [ 2*T_b_s1e*Pr_st1X  zeros(3,12)  Pc1X
           2*T_b_c*PrX        zeros(3,12)  zeros(3,3) ];        % sensitivity matrix % VECTOR    
    %H  = [ T_b_s1e           zeros(3,12) .5*eye(3)     % doesn't work
    %       2*[cross(T_b_c(1,:))*Pr]' zeros(1,12) zeros(1,3) 
    %       2*[cross(T_b_c(1,:))*Pr]' zeros(1,12) zeros(1,3)
    %       2*[cross(T_b_c(1,:))*Pr]' zeros(1,12) zeros(1,3) ];        % sensitivity matrix
    %H  = [ T_b_s1e     zeros(3,12) .5*eye(3) ];         % without payload
    
    % update
    K = P*H'*inv(H*P*H'+R);
    P = (eye(18)-K*H)*P;
    dx = K*z; dq=dx(1:3); db=dx(4:6); dg=dx(7:15); ds1=dx(16:18);
    ge=ge+dg; 
    %we = T_g0_b*[(wgm+be) + (eye(3)+Me)*db + omega_g(wgm+be)*ge];      % to avoid linearization at early stage
    %qe=qe+xi(qe)*dq; qe=qe/norm(qe);
    dq=[dx(1:3);1];  qe=q_mult(dq,qe); qe=qe/norm(qe);      
    be=be+db; s1e=s1e+ds1;
    %we = T_g0_b*[wgm + (eye(3)+Me)*be + omega_g(wgm)*ge];
    qe_his(:,i+1)=qe; we_his(:,i+1)=we; be_his(:,i+1)=be; ge_his(:,i+1)=ge; 
    s1e_his(:,i+1)=s1e; P_cov(:,i+1)=diag(P);
    %qe=expm(.5*dts*omega(we))*qe;       %%%%% !!!!!  % this one works good    
    waitbar(i/(tf/dts),wait)
end
close(wait)

% these are the parameters to be sent back
qe=qe_his; we=we_his; be=be_his; ge=ge_his; s1e=s1e_his;

